use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::Write;
use rand_distr::WeightedAliasIndex;

use rayon::slice::ParallelSliceMut;

use crate::cedar::Cedar;
use crate::cedar::readers::{Query, Mapping};
use rand::prelude::*;

mod get_taxonomy;
use get_taxonomy::tax_main;

/*
    Candidate move for an incumbent to make room on a reference
 */
struct CandidateMove {
    cost: f32,
    incumbent_id: usize,
    incumbent_score: usize,
    alt_ref_id: usize,
    alt_score: usize
}

/*
    abundance: <ref_id, estimated abundance from cedar>
    current_abundance: <ref_id, current abundance>
    assignments: <ref_id, <score, query_name>>
    output_assignments: <query_name, ref_id>
    query_size: the total number of queries (not counting unmapped queries)
    max_diff: maximum allowed difference between abundance and current abundance
*/
struct AssignmentMachine {
    abundance: HashMap<usize, f32>,
    current_abundance: HashMap<usize, f32>,
    assignments: HashMap<usize, HashMap<usize, HashSet<usize>>>,
    output_assignments: HashMap<usize, usize>,
    query_size: usize,
    max_diff: f32,
}

impl AssignmentMachine {
    fn new(abundance: HashMap<usize, f32>, query_size: usize) -> AssignmentMachine {
        let mut current_abundance = HashMap::new();
        for (key, _) in &abundance {
            current_abundance.insert(*key, 0.0);
        }
        AssignmentMachine { abundance, current_abundance, assignments: HashMap::new(), output_assignments: HashMap::new(), query_size, max_diff: 1.0/(query_size as f32)}
    }

    // assign a query to a reference with a score
    fn add_assignment(&mut self, query_id: usize, ref_id: usize, score: usize) {
        if ref_id != usize::MAX - 1 {
            let entry = self.current_abundance.entry(ref_id).or_insert(0.0);
            *entry += 1.0 / self.query_size as f32;                             // update the new current abundance
        }
        self.output_assignments.insert(query_id, ref_id);

        let entry = self.assignments.entry(ref_id).or_insert(HashMap::new());
        let entry2 = entry.entry(score).or_insert(HashSet::new());
        entry2.insert(query_id);
    }

    // check if a reference still has space to be assigned a query
    fn has_space(&self, ref_id: &usize) -> bool {
        if *ref_id == usize::MAX - 1 {
            true
        } else {
            self.current_abundance[ref_id] < self.abundance[ref_id] + self.max_diff
        }
    }

    // remove an assignment for a previously assigned query
    fn remove_assignment(&mut self, ref_id: usize, query_id: usize, score: usize) {
        let entry = self.current_abundance.entry(ref_id).or_insert(0.0);
        *entry -= 1.0 / self.query_size as f32;
        
        let entry = self.assignments.entry(ref_id).or_insert(HashMap::new());
        let entry2 = entry.entry(score).or_insert(HashSet::new());
        entry2.remove(&query_id);
    }

    // initial assignment of unique mapping queries or unmapped queries
    fn initial_assignment(&mut self, mut queries: HashMap<usize, Query>) -> HashMap<usize, Query> {
        let mut unique_mapping_queries = 0;
        let mut unmapped_queries = 0;
        let mut to_remove = Vec::new();
        
        println!("total query_size = {}", self.query_size);
        for (query_id, query) in &queries {
            if query.mappings.is_empty() {
                unmapped_queries += 1;

                self.add_assignment(*query_id, usize::MAX, 0);
                to_remove.push(*query_id)
            } else if check_if_read_is_uniquely_mapping(&query.mappings) {
                unique_mapping_queries += 1;

                let ref_id = query.mappings.iter().next().unwrap().get_reference_id();
                self.add_assignment(*query_id, ref_id, 1000);
                to_remove.push(*query_id)
            }
        }

        for id in to_remove {
            queries.remove(&id);
        }

        println!("unique mapping queries: {}, unmapped queries: {}, current assigned length: {}", unique_mapping_queries, unmapped_queries, self.output_assignments.len());

        queries
    }

    // secondary assignment of queries whose best mapping is a lot better than the secondary mappings
    fn secondary_assignment(&mut self, mut queries: HashMap<usize, Query>, score_max_diff: f32) -> HashMap<usize, Query>{
        let mut to_remove = Vec::new();

        for (query_id, query) in &queries {
            let (ref_id, score) = assign_best(query, &self, score_max_diff);
            if ref_id != usize::MAX - 1 {
                to_remove.push(*query_id);
                if self.has_space(&ref_id) {
                    self.add_assignment(*query_id, ref_id, score as usize);
                }
            }
        }

        for id in to_remove {
            queries.remove(&id);
        }

        queries
    }

    fn compute_candidate_moves(
        &self,
        original_queries: &HashMap<usize, Query>,
        leftovers: &HashMap<usize, Query>
    ) -> HashMap<usize, VecDeque<CandidateMove>>
    {
        let mut max_thresholds: HashMap<usize, f32> = HashMap::new();
        for query in leftovers.values() {
            let total = query.get_total_score();

            for mapping in &query.mappings {
                let new_threshold = mapping.get_score() / total;
                let cur_threshold = max_thresholds
                    .entry(mapping.get_reference_id())
                    .or_insert(f32::NEG_INFINITY);

                if new_threshold > *cur_threshold {
                    *cur_threshold = new_threshold;
                }
            }
        }

        let mut candidates = HashMap::new();
        for (&ref_id, &threshold_limit) in &max_thresholds {
            let bins = match self.assignments.get(&ref_id) {
                Some(b) => b,
                None => continue
            };
            let mut list = Vec::new();

            for (score, incumbents) in bins {
                for inc_id in incumbents {
                    let inc = &original_queries[inc_id];
                    let total_score = inc.get_total_score();

                    for mapping in &inc.mappings {
                        let alt_ref_id = mapping.get_reference_id();

                        if alt_ref_id != ref_id && self.has_space(&alt_ref_id) {
                            let cost = 2.0 * (*score as f32) / total_score
                                - mapping.get_score() / total_score;

                            if cost < threshold_limit {
                                list.push(CandidateMove {
                                    cost,
                                    incumbent_id: *inc_id,
                                    incumbent_score: *score,
                                    alt_ref_id,
                                    alt_score: mapping.get_score() as usize
                                })
                            }
                        }
                    }
                }
            }

            if !list.is_empty() {
                list.sort_by(|a, b| a.cost.partial_cmp(&b.cost).unwrap());
                candidates.insert(ref_id, VecDeque::from(list));
            }
        }

        candidates
    }

    fn place_leftovers(&mut self, leftovers: &HashMap<usize, Query>,
        candidates: &mut HashMap<usize, VecDeque<CandidateMove>>)
        -> HashMap<usize, Query>
    {
        let mut unplaced = HashMap::new();

        let mut leftover_ids: Vec<usize> = leftovers.keys().cloned().collect();
        leftover_ids.sort_by(|a, b| {
            leftovers[b].get_best_mappings().1
                .partial_cmp(&leftovers[a].get_best_mappings().1).unwrap()
        });

        for query_id in leftover_ids {
            let total_score = leftovers[&query_id].get_total_score();
            let mut placed = false;

            for mapping in leftovers[&query_id].sort_mappings() {
                let target = mapping.get_reference_id();
                let threshold = mapping.get_score() / total_score;

                let chosen = match candidates.get_mut(&target) {
                    Some(queue) => self.take_valid_move(queue, target, threshold),
                    None => continue
                };

                if let Some(c) = chosen {
                    self.remove_assignment(target, c.incumbent_id, c.incumbent_score);
                    self.add_assignment(c.incumbent_id, c.alt_ref_id, c.alt_score);
                    self.add_assignment(query_id, target, mapping.get_score() as usize);
                    placed = true;
                    break;
                }
            }

            if !placed {
                unplaced.insert(query_id, leftovers[&query_id].clone());
            }
        }

        unplaced
    }

    fn take_valid_move(&self, queue: &mut VecDeque<CandidateMove>, target: usize, threshold: f32)
        -> Option<CandidateMove>
    {
        while !queue.is_empty() {
            if queue[0].cost >= threshold {
                return None;
            }

            let (inc_id, inc_score, alt_ref_id) = (
                queue[0].incumbent_id, queue[0].incumbent_score, queue[0].alt_ref_id
            );

            let incumbent_here = self.assignments.get(&target)
                .and_then(|b| b.get(&inc_score))
                .is_some_and(|s| s.contains(&inc_id));

            if incumbent_here && self.has_space(&alt_ref_id) {
                return queue.pop_front()
            }
            queue.pop_front();
        }

        None
    }

    //assignment based on abundances
    fn assign_based_on_abundance(&mut self, mut queries: HashMap<usize, Query>, original_queries: &HashMap<usize, Query>, method: String) {
        //create score bins
        let mut score_bins = HashMap::new();
        for (query_id, query) in &queries {
            for mapping in &query.mappings {
                if self.has_space(&mapping.get_reference_id()) {
                    let entry = score_bins.entry(mapping.get_score() as usize).or_insert(Vec::new());
                    entry.push((*query_id, mapping.get_reference_id()));
                }
            }
        }
        //sort the bins from high scores to low scores
        let mut keys: Vec<usize> = score_bins.keys().cloned().collect();
        keys.par_sort_by(|a, b| b.cmp(a));
        
        // map the queries in order
        for key in keys {
            for (query_id, ref_id) in &score_bins[&key] {
                if queries.contains_key(query_id) && self.has_space(ref_id){
                    self.add_assignment(*query_id, *ref_id, key);
                    queries.remove(query_id);
                }
            }
        }
        println!("queries have been mapped to the best possible reference. left overs: {}", queries.len());

        // assignment of the left overs by trying to open up space
        let leftover_queries = if queries.is_empty() {
            HashMap::new()
        } else {
            let mut candidates = self.compute_candidate_moves(original_queries, &queries);
            self.place_leftovers(&queries, &mut candidates)
        };
        
        println!("cannot move: {}", leftover_queries.len());

        //assign the queries that cannot be mapped based on what method was specified. Default mode is none
        if method == "none" {
            self.leave_left_overs(leftover_queries);
        } else {
            self.assign_based_on_prob(leftover_queries);
        }
    }

    //assignment based on probability with weights being the mapping scores
    fn assign_based_on_prob(&mut self, queries: HashMap<usize, Query>) {
        for (name, query) in queries {
            let mappings:Vec<Mapping> = query.mappings.into_iter().collect();
            let dist  = WeightedAliasIndex::new(mappings.iter().map(|mapping| mapping.get_score()).collect()).unwrap();
            // assign randomly
            let mut rng = thread_rng();
            let chosen = &mappings[dist.sample(&mut rng)];
            self.add_assignment(name, chosen.get_reference_id(), chosen.get_score() as usize);
        }
    }

    // leave left over queries unassigned
    fn leave_left_overs(&mut self, queries: HashMap<usize, Query>) {
        for (name, _query) in queries {
            self.add_assignment(name, usize::MAX - 1, 0);
        }
    }
}


// helper to check if read is uniquely mapping
fn check_if_read_is_uniquely_mapping(mappings: &HashSet<Mapping>) -> bool {
    let mut it = mappings.iter();
    let first = match it.next() {
        Some(m) => m.get_reference_id(),
        None => return false
    };
    it.all(|m| m.get_reference_id() == first)
}

// assign each mapping to a unique reference based on their mapping scores and the predicted abundance levels
pub(crate) fn assign_mappings(cedar: Cedar, score_max_diff: f32, method: String) -> HashMap<String, String> {
    let mut queries = cedar.get_queries();
    let references = cedar.get_references();
    let mut machine = AssignmentMachine::new(cedar.get_strain_abundance(), queries.len() - cedar.get_unmapping_reads());

    println!("performing assignment of queries\n");

    //initial assignment
    queries = machine.initial_assignment(queries);

    println!("initial assignment done. query length: {}", queries.len());
    
    //secondary assignment
    queries = machine.secondary_assignment(queries, score_max_diff);

    println!("second assignment done. query length: {}", queries.len());

    //assignment of the rest of the queries based on abundancies
    machine.assign_based_on_abundance(queries, &cedar.queries, method);

    println!("final assignment done. output length: {}", machine.output_assignments.len());

    //write final assignments
    let mut output = HashMap::new();
    for (query_id, id) in machine.output_assignments {
        let name = cedar.query_id_2_name[&query_id].to_string();
        if id == usize::MAX - 1 {                                           // queries that could not be assigned to anything (if using method 1: leave left over reads un-assigned)
            output.insert(name, "NOT ALIGNED".to_string());
        } else if id == usize::MAX {                                        // queries that couldn't be mapped to any reference from the first aligner
            output.insert(name, "NOT ALIGNED".to_string());  
        } else {                                                            // queries that were assigned to something
            output.insert(name, references[&id].ref_name.to_string());
        }
    }
    output
}

// find the mapping with the best score and if it is a lot bigger than the second best mapping, return it
fn assign_best(query: &Query, machine: &AssignmentMachine, score_max_diff: f32) -> (usize, f32) {
    let (r1, s1, _r2, s2) = query.get_best_mappings();
    if (s1 - s2)/s1 > score_max_diff {
    // difference is enough to make it get mapped to the best mapping
        if  machine.has_space(&r1){
            (r1, s1)
        } else {
            (usize::MAX - 1, 0.0)
        }
    } else {
        (usize::MAX - 1, 0.0)
    }
}

// write the output into a file in the following way: query_name    reference_name
pub(crate) fn write_output(output_filename: String, output: HashMap<String, String>) {
    let mut output_file = File::create(output_filename).unwrap();
    for (q_name, r_name) in output {
        let mut data = q_name;
        data.push_str("\t");
        data.push_str(&r_name);
        data.push_str("\t"); 
        data.push_str("\n");
        output_file.write(data.as_bytes()).ok();
    }
}

// write the output into a file in the following way: query_name    reference_name  reference_species   reference_genus     reference_family    ...     reference_superkingdom
pub(crate) fn write_output_with_taxonomy(output_filename: String, output: HashMap<String, String>, at_file: String, nodes_file: String, names_file: String) {
    tax_main(output, at_file, nodes_file, names_file, output_filename)
}
