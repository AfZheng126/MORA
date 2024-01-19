use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::sync::Mutex;
use std::sync::mpsc::channel;
use rand_distr::WeightedAliasIndex;

use rayon::prelude::{IntoParallelIterator, ParallelIterator, IntoParallelRefIterator};
use rayon::slice::ParallelSliceMut;

use crate::cedar::Cedar;
use crate::cedar::readers::{Query, Mapping};
use rand::prelude::*;

mod get_taxonomy;
use get_taxonomy::tax_main;
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
    fn initial_assignment(&mut self, queries: HashMap<usize, Query>) -> HashMap<usize, Query> {
        let mut updated_queries = queries.clone(); 

        let mut unique_mapping_queries = 0;
        let mut unmapped_queries = 0;
        
        println!("total query_size = {}", self.query_size);
        for (query_id, query) in queries {
            if query.mappings.len() == 1 {
                unique_mapping_queries += 1;

                updated_queries.remove(&query_id);
                let maps: Vec<Mapping> = query.mappings.into_iter().collect();
                self.add_assignment(query_id, maps[0].get_reference_id(), 1000);
            } else if query.mappings.len() == 0 {
                unmapped_queries += 1;
                updated_queries.remove(&query_id);
                self.add_assignment(query_id, usize::MAX, 0);
            }
        }
        println!("unique mapping queries: {}, unmapped queries: {}, current assigned length: {}", unique_mapping_queries, unmapped_queries, self.output_assignments.len());
        updated_queries
    }

    // secondary assignment of queries whose best mapping is a lot better than the secondary mappings
    fn secondary_assignment(&mut self, queries: HashMap<usize, Query>, score_max_diff: f32) -> HashMap<usize, Query>{
        let mut updated_queries = queries.clone();
        for (query_id, query) in &queries {
            let (ref_id, score) = assign_best(query, &self, score_max_diff);
            if ref_id != usize::MAX - 1 {
                updated_queries.remove(query_id);
                if self.has_space(&ref_id) {
                    self.add_assignment(*query_id, ref_id, score as usize);
                }
            }
        }
        updated_queries
    }

    //assignment based on abudancies
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
        
        // assignment of the left overs by tring to open up space
        let cannot_move = Mutex::new(0);
        let left_over_queries = Mutex::new(HashMap::new());
        if queries.len() != 0 {
            let (sender, receiver) = channel();

            queries.par_iter().for_each_with(sender, |s, (query_id, query)|{
                let ordered_mappings = query.sort_mappings();
                let total_score = query.get_total_score();
                //check if something can be moved from the reference that the mapping is mapped to
                let mut moved = false;
                for map in ordered_mappings {
                    let new_score = map.get_score();
                    let new_normalized_score = map.get_score() / total_score;
                    let (q_moved, id_of_moved, re_assigned_ref, original_score, re_assigned_score) = try_open_up_space(&self.assignments[&map.get_reference_id()], self, &map.get_reference_id(), original_queries, new_normalized_score);
                    if q_moved {
                        s.send((map.get_reference_id(), id_of_moved, original_score, re_assigned_ref, re_assigned_score, query_id, map.get_reference_id(), new_score as usize)).ok();
                        moved = true;
                        break;
                    }
                }
                if !moved {
                    let mut locked = cannot_move.lock().unwrap();
                    *locked += 1;
                    let mut left_over_queries_locked = left_over_queries.lock().unwrap();
                    left_over_queries_locked.insert(*query_id, query.clone());
                }
            });

            let mut res: Vec<_> = receiver.iter().collect();
            res.iter_mut().for_each(|x| {
                let (old_ref_id, id_of_moved, original_score, new_ref_id, new_score, 
                    query_id, map_ref_id, map_score) = x;
                self.remove_assignment(*old_ref_id, *id_of_moved, *original_score);
                self.add_assignment(*id_of_moved, *new_ref_id, *new_score);
                self.add_assignment(**query_id, *map_ref_id, *map_score);
            });
        }
        println!("cannot move: {}", cannot_move.into_inner().unwrap());

        //assign the queries that cannot be mapped based on what method was specified. Default mode is none
        if method == "none" {
            self.leave_left_overs(left_over_queries.into_inner().unwrap());
        } else {
            self.assign_based_on_prob(left_over_queries.into_inner().unwrap());
        }
    }

    //assignment based on probability with weights being the mapping scores
    fn assign_based_on_prob(&mut self, queries: HashMap<usize, Query>) {
        for (name, query) in queries {
            let mappings:Vec<Mapping> = query.mappings.into_par_iter().collect();
            let dist  = WeightedAliasIndex::new(mappings.par_iter().map(|mapping| mapping.get_score()).collect()).unwrap();
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

// assign each mapping to a unique reference based on their mapping scores and the predicted abundance levels
pub(crate) fn assign_mappings(cedar: Cedar, score_max_diff: f32, method: String) -> HashMap<String, String> {
    let mut queries = cedar.get_queries();
    let references = cedar.get_references();
    let mut machine = AssignmentMachine::new(cedar.get_strain_abundance(), queries.len() - cedar.get_unmapping_reads());

    println!("performing assignment of quries\n");

    //initial assignment
    queries = machine.initial_assignment(queries);

    println!("intial assignment done. query length: {}", queries.len());
    
    //secondary assignment
    queries = machine.secondary_assignment(queries, score_max_diff);

    println!("second assignment done. query length: {}", queries.len());

    //assignment of the rest of the queries based on abundancies
    machine.assign_based_on_abundance(queries, &cedar.get_queries(), method);

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

/*
    assignemnts: <score, query_name>
    ref_id: the reference that we are trying to open up space in 
    new_normalized_score: the score for the read we are trying to make space for

    output: 
        - if it has been moved
        - id of query that has been moved
        - id of reference query was moved to
        - original mapping score
        - new mapping score
*/
fn try_open_up_space(assignments: &HashMap<usize, HashSet<usize>>, machine: &AssignmentMachine, ref_id: &usize, queries: &HashMap<usize, Query>, new_normalized_score: f32) -> (bool, usize, usize, usize, usize) {
    if assignments.len() == 0 {
        return (false, 0, 0, 0, 0);
    }

    // sort the keys from smallest to largest
    let mut keys: Vec<&usize> = assignments.keys().collect();
    keys.par_sort_by(|a, b| a.cmp(b));

    // try find a query that can be moved
    for key in &keys {
        for query_id in assignments.get(key).unwrap() {
            let query = &queries[query_id];
            let total_score = query.get_total_score();
            // try assign the query to somewhere else
            for mapping in &query.mappings {
                let re_assigned_ref_id = &mapping.get_reference_id();
                if re_assigned_ref_id != ref_id && machine.has_space(re_assigned_ref_id) {
                    if (**key as f32)/total_score - mapping.get_score()/total_score  <  new_normalized_score - (**key as f32)/total_score {
                        return (true, *query_id, *re_assigned_ref_id, **key, mapping.get_score() as usize);
                    }
                }
            }
        }
    }
    (false, 0, 0, 0, 0)
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
