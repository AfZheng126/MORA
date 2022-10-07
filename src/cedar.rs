use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;

mod equivalence_class_builder;
use equivalence_class_builder::EquivalenceClassBuilder;

mod stats;
use stats::Stats;

mod taxa;
use taxa::TaxaNode;

mod set_covers;
use set_covers::greedy_set_cover;

pub(crate) mod readers;
use readers::{read_by_bash, Query, Reference};

use crate::cedar::equivalence_class_builder::TGValue;
use crate::cedar::readers::{read_rank, read_taxonomy};

use self::equivalence_class_builder::TargetGroup;

mod util;

use std::sync::Mutex;
use rayon::prelude::*;
use std::sync::mpsc::channel;

// helper functions 
fn create_eqb(read_per_strain_prob_inst: Vec<(usize, f32)>) -> (TargetGroup, Vec<f32>) {
    // sort the read_per_strain_prob_inst so that the usize is in order of smallest to biggest
    let sorted_vec = util::sort_vec(&read_per_strain_prob_inst);
    let mut genome_ids = Vec::new();
    genome_ids.reserve(read_per_strain_prob_inst.len());
    let mut probs = Vec::new();
    probs.reserve(read_per_strain_prob_inst.len());

    for it in sorted_vec {
        genome_ids.push(it.0);
        probs.push(it.1);
    }
    let tgt = TargetGroup::new_with_tgt(genome_ids);
    (tgt, probs)
}


#[derive(Clone)]
pub(crate) struct Cedar {
    /*
    eqb: The equivalence Class Builder
    taxa_node_map: 
    ref_name_2_tax_id: contains the names of the references and their taxonomy level
    ref_id_2_tax_id: contains the ids of the reference and their taxonomy level 
    re_id_2_name: contains the ids of the reference and their name in the SAM file
    query_id_2_name: contains the ids of hte queries and their name in the SAM file

    strain_coverage: coverage of each reference
    strain_coverage_bins: the references IDs and their corresponding bin coverages 
    strain_abundance: contains the ref_IDs and their percent abundance
    taxa_abundance: contains the taxa_IDs and their abundance
    unampping_reads: the number of queries that are mapped to 0 references

    read_cnt: the total number of reads
    cov: the coverage of each reference 
    queries: the hashmap of query names to queries
    references: the hashmap of ref_ids to references, ref_id starts at 1
    */
    eqb: EquivalenceClassBuilder,
    taxa_node_map: HashMap<usize, TaxaNode>,
    ref_name_2_tax_id: HashMap<String, usize>,
    ref_id_2_name: HashMap<usize, String>,
    pub(crate) query_id_2_name: HashMap<usize, String>,
    strain_coverage: HashMap<usize, f32>,
    strain_coverage_bins: HashMap<usize, Vec<usize>>, 
    strain_abundance: HashMap<usize, f32>, 
    read_cnt: usize,
    ref_id_to_tax_id: HashMap<usize, usize>,
    cov: HashMap<usize, usize>,
    pub(crate) queries: HashMap<usize, Query>,
    references: HashMap<usize, Reference>,
    taxa_abundance: HashMap<usize, f32>,
    unmapping_reads: usize,
}

impl Cedar {
    pub(crate) fn new(taxonomy_tree_filename: String, ref_name_2_tax_id_filename: String, 
        flat_abundance: bool, batch_size: usize) -> Cedar {
        println!{"Constructing Cedar"};

        let mut ref_name_2_tax_id = HashMap::new();
        let mut taxa_node_map = HashMap::new();
        if !flat_abundance {
            // map rank string values to enum values
            ref_name_2_tax_id = read_rank(ref_name_2_tax_id_filename, batch_size);
            taxa_node_map = read_taxonomy(taxonomy_tree_filename, batch_size);
        }
        
        Cedar { eqb: EquivalenceClassBuilder::new(), taxa_node_map, ref_name_2_tax_id, strain_coverage: HashMap::new(), 
            strain_coverage_bins: HashMap::new(), strain_abundance: HashMap::new(), read_cnt: 0, ref_id_to_tax_id: HashMap::new(), ref_id_2_name: HashMap::new(), query_id_2_name: HashMap::new(), cov: HashMap::new(), 
            queries: HashMap::new(), references: HashMap::new(), taxa_abundance: HashMap::new(), unmapping_reads: 0 }
    }

    // find the stats of the current list of queries and also updates the equivalence class builder
    fn process_reads_parallel(&mut self, flat_abundance: bool) -> Stats {
        let (total_read_cnt, total_multi_mapped_reads, total_unmapped_reads) = (Mutex::new(0), Mutex::new(0), Mutex::new(0));
        let (sender, receiver) = channel();
 
        self.queries.par_iter().for_each_with(sender, |s, (_name, query)| {
            if query.get_cnt() != 0 {
                if query.get_cnt() > 1 {
                    *total_multi_mapped_reads.lock().unwrap() += 1;
                }
                let mut read_per_strain_prob_inst = Vec::new(); 
                let mut mapping_score = 0;
                for mapping in &query.mappings {
                    mapping_score += mapping.get_score() as usize;
                    let value = mapping.get_score() / self.references.get(&mapping.get_reference_id()).unwrap().ref_len as f32;
                    read_per_strain_prob_inst.push((mapping.get_reference_id(), value));
                }
                s.send((mapping_score, read_per_strain_prob_inst)).ok();
            } else {
                *total_unmapped_reads.lock().unwrap() += 1;
            }
            *total_read_cnt.lock().unwrap() += 1;
        });

        let mut res: Vec<_> = receiver.iter().collect();

        res.iter_mut().for_each(|x| {
            let mapping_score = x.0;
            let read_per_strain_prob_inst = &x.1;
            for i in 0..read_per_strain_prob_inst.len() {
                let entry = self.strain_abundance.entry(read_per_strain_prob_inst[i].0).or_insert(0.0);
                *entry += 1.0 / read_per_strain_prob_inst.len() as f32;

                let tid;
                if flat_abundance {
                    tid = read_per_strain_prob_inst[i].0;
                } else {
                    let ref_name = &self.references[&read_per_strain_prob_inst[i].0].ref_name;
                    tid = self.ref_name_2_tax_id[ref_name];
                }
                let entry = self.cov.entry(tid).or_insert(0);
                *entry += mapping_score;
            }

            let c = create_eqb((*x.1).to_vec());
            self.eqb.add_group(c.0, c.1);
        });

        Stats::new_with_stats(total_read_cnt.into_inner().unwrap(), total_multi_mapped_reads.into_inner().unwrap(), total_unmapped_reads.into_inner().unwrap())
    }


    // updates the bins 
    // the mappings in the self.queries[i] are the ones that have a reference_id, so we don't need to check if the reference_id is valid
    fn update_bins(&mut self, segment_size: usize) {
        for (_, query) in &self.queries {
            for mapping in &query.mappings {
                let bin_number = mapping.get_position() as usize / segment_size;
                if mapping.get_reference_id() == usize::MAX {
                    continue;
                }
                let entry = self.strain_coverage_bins.entry(mapping.get_reference_id()).or_insert(Vec::new());
                entry[bin_number as usize] += 1;
            }
        }
    }

    // updates the coverage of the references 
    fn calculate_coverage(&mut self) {
        for i in self.strain_coverage_bins.keys() {
            let bins = self.strain_coverage_bins.get(&i).unwrap();
            let mut covered = 0.0;
            // let mut expression_depth = 0;
            for j in 0..bins.len() {
                if bins[j] > 0 {
                    covered += 1.0;
                }
                // expression_depth += bins[j];
            }
            self.strain_coverage.insert(*i, covered / bins.len() as f32 );
        }
    }

    /*
    read a SAM file and extract the information from into and store it in the desired form in the Cedar struct    
    Inputs:
    mapper_output_filename: The path to the SAM file
    require_concordance:
    flat_abundance: 
    only_unique:
    only_perfect:
    segment_size: the size of each bin
    range_factorization: 
    */
    pub(crate) fn load_mapping_info_parallel(&mut self, mapper_output_filename: String, flat_abundance: bool, segment_size: usize, batch_size: usize, method: String) {
        println!("Cedar: Load Mapping File");
        println!("Mapping Ouput File: {}", mapper_output_filename);
        
        // load the information from the file
        let c = read_by_bash(mapper_output_filename, batch_size, method);
        self.references = c.0;
        self.queries = c.1;
        self.query_id_2_name = c.2;

        for key in self.references.keys() {
            self.strain_abundance.insert(*key, 0.0);
        }

        // constucting coverage bins;
        let (sender, receiver) = channel();

        (0..self.references.len()).into_par_iter().for_each_with(sender, |s, index| {
            let ref_name = &self.references.get(&(index + 1)).unwrap().ref_name;
            let ref_length = self.references.get(&(index + 1)).unwrap().ref_len;
            let bin_cnt = ref_length / segment_size + 1;
            let bins = vec![0; bin_cnt];
            let tid;
            if flat_abundance {
                tid = index + 1;
            } else {
                tid = self.ref_name_2_tax_id.get(ref_name).unwrap().to_owned();
            }
            let data = (index + 1, bins, tid);
            s.send(data).ok();
        });

        let res: Vec<_> = receiver.iter().collect();
        
        res.iter().for_each(|x| {
            let index = x.0;
            let bins = &x.1;
            let tid = x.2;
            self.strain_coverage_bins.insert(index, bins.to_vec());
            self.ref_id_to_tax_id.insert(index, tid);
            self.cov.insert(tid, 0);
        });

        // update the information in the Cedar struct
        let stats = self.process_reads_parallel(flat_abundance);
        self.read_cnt = stats.get_total_read_cnt();
        self.update_bins(segment_size);
        self.calculate_coverage();
        self.unmapping_reads = stats.get_total_unmapped_reads();

        // print the information obtained from the mapping_output_file 
        stats.print_stats();
    }

    // applying the (greedy) set cover to find the minimum number of references that covers all the equivalence classes
    fn apply_set_cover(&self, strain_cnt: &Vec<f32>, mut strain_valid: HashMap<usize, bool>, 
        mut strain_potentially_removable: HashMap<usize, bool>, min_cnt: f32, mut can_help: bool) -> (bool, HashMap<usize, bool>, HashMap<usize, bool>) {

        let mut previously_valid:u64 = 0;
        let eq_vec = &self.eqb.count_vec;

        // rule to find potentially removable strains
        for i in 1..strain_cnt.len() {
            if strain_cnt[i] <= min_cnt {
                strain_potentially_removable.insert(i, strain_valid[&i]);
            } else {
                strain_potentially_removable.insert(i, false);
            }
            if strain_valid[&i] == true {
                previously_valid += 1;
            }
        }
        // find the list of all references with a unique "valid" read mapped to them
        let mut unique_reads_refs = HashSet::new();
        for i in 0..eq_vec.len() {
            let tg = &eq_vec[i].0;
            let v = &eq_vec[i].1;
            let csize = v.get_weights().len();
            let mut eq_valids = Vec::new();
                
            for read_mapping_cntr in 0..csize {
                let tgt = tg.get_tgts()[read_mapping_cntr]; 
                if strain_valid[&tgt] {
                    eq_valids.push(tgt);
                }
            }
            if eq_valids.len() == 1 {
                unique_reads_refs.insert(eq_valids[0]);
            }
        }
        // if the reference has only one query mapping to it, it is potentially removable
        for i in 1..strain_cnt.len() {
            if strain_valid[&i] && unique_reads_refs.contains(&i) != true {
                strain_potentially_removable.insert(i, true);
            }
        }
        
        // ref_2_eqset is a map from a reference to the set of ambiguous eqs it belongs to.
        // ambiguous eq is an equivalence class such that all its remaining refs has been chosen as potentially removable. 
        let mut ref_2_eqset:HashMap<usize, HashSet<u64>> = HashMap::new();       // key: reference_ID, value: set of ambiguous equivalence classes
        
        for i in 0..eq_vec.len() {                
            let tg = eq_vec[i].0.clone();
            let v =  eq_vec[i].1.clone();
            let csize = v.get_weights().len();
            let mut total_valid_refs_in_cur_eq = 0;
            let mut potential_removable_cntr = 0;
    
            // count total number of valid reference for each eq.
            // (other than those that have been invalidated in previous rounds)
            for read_mapping_cntr in 0..csize {
                let tgt = tg.get_tgts()[read_mapping_cntr];
                if strain_valid[&tgt] {
                    total_valid_refs_in_cur_eq += 1;
                }
                if strain_potentially_removable[&tgt] {
                    potential_removable_cntr += 1;
                }
            }

            // if all the refs in the eq are set as potentially_removable then it is an ambiguous eq and we add it to ref2_eqset
            if potential_removable_cntr >= total_valid_refs_in_cur_eq { 
                for read_mapping_cntr in 0..csize {
                    let tgt = tg.get_tgts()[read_mapping_cntr];
                    if strain_valid[&tgt] && strain_potentially_removable[&tgt] {
                        let entry = ref_2_eqset.entry(tgt).or_insert(HashSet::new());
                        entry.insert(tg.get_hash());
                    }
                }
            }
        }
    
        let mut eq_2_id:HashMap<u64, u64> = HashMap::new();         // key: reference ID, value: an ID we use for the greedy algorithm

        if ref_2_eqset.len() > 0 {
            // set cover input preparation
            // convert the input to proper format for the library that runs setCover algo. 
            let mut id = 1;
            // turns all the references in ref_2_eqset into a list
            for kv in ref_2_eqset.keys() {
                for v in &ref_2_eqset[kv] {
                    if eq_2_id.contains_key(&v) != true { 
                        eq_2_id.insert(*v, id);
                        id += 1;
                    }
                }
            }

            let mut sets = Vec::new();
            let mut weights = Vec::new();
            let unique_element_count = eq_2_id.len();
            let mut set_cover_id_2_ref = Vec::new();
            let mut set_weight = Vec::new();

            for kv in ref_2_eqset {
                set_cover_id_2_ref.push(kv.0);
                let set_size = kv.1.len();
                let mut mw = (self.strain_coverage[&kv.0] * 10000.0) as usize;
                let mut per_element = mw / set_size;
                if per_element < 1 {
                    per_element = 1;
                }
                mw = per_element * set_size;
                set_weight.push(mw);

                let mut set = Vec::new();
                for eq in kv.1 {
                    set.push(eq_2_id[&eq]);
                }
                sets.push(set);
                weights.push(vec![per_element as f32; set_size]);
            }

            // end of set cover input preparation
            // run set_cover algorithm over the lists of refs and eqs in ref2_eqset
            // as we are givng sets and weights and not &sets and &weights, they are lost after greedy_set_cover finishes
            let final_covering = greedy_set_cover(sets, weights, unique_element_count);

            // put the list of minimum # of references that can cover all eqs in remainingRefs
            let mut remaining_refs = HashSet::new();
            for it in final_covering {
                remaining_refs.insert(set_cover_id_2_ref[it]);
            }

            // go over the list of references
            for ref_cntr in 1..(strain_valid.len() + 1) {
                if strain_potentially_removable[&ref_cntr] && remaining_refs.contains(&ref_cntr) != true {
                    strain_valid.insert(ref_cntr, false);
                }
            }

            let mut total_valid = 0;
            for s in 0..strain_valid.len() {
                if strain_valid[&(s+1)] {
                    total_valid += 1;
                }
            }
            can_help = previously_valid != total_valid;
        }
        (can_help, strain_valid, strain_potentially_removable)
    }

    // the EM function to calculate the abundances
    pub(crate) fn parallel_em(&mut self, max_iter: usize, eps: f32, min_cnt: f32) {
        self.eqb.finish();

        let eq_vec:Vec<(TargetGroup, TGValue)> = self.eqb.count_vec.par_iter().cloned().collect();

        // finds the maximum sequence ID in the strains.
        let max_seq_id = self.strain_abundance.len();

        let mut new_strain_cnt = vec![0.0; max_seq_id + 1];
        //let mut strain_cnt:Vec<f32> = vec![0.0; (max_seq_id + 1) as usize];
        let mut strain_valid:HashMap<usize, bool> = HashMap::new();
        let mut strain_potentially_removable:HashMap<usize, bool> = HashMap::new();

        for i in 1..(max_seq_id + 1) {
            strain_valid.insert(i, true);
            strain_potentially_removable.insert(i, false);
        }

        //clone the data in self.strain
        let mut strain_cnt = vec![0.0; max_seq_id + 1];
        for (key, value) in &self.strain_abundance {
            strain_cnt[*key] = *value;
        }

        let tot_count: usize = eq_vec.par_iter().map(|x| x.1.get_count()).sum();

        // logger stuff
        println!("Total starting count {}", tot_count);
        println!("Total mapped reads cnt {}", self.read_cnt); 

        let mut cntr:usize = 0;
        let mut converged = false;
        let thresholding_iter_step = 10;
        let mut can_help = true;

        while cntr < max_iter && converged == false {
            if cntr % thresholding_iter_step == 0 && can_help {
                let a = self.apply_set_cover(&strain_cnt, strain_valid, strain_potentially_removable, min_cnt, can_help);
                can_help = a.0;
                strain_valid = a.1;
                strain_potentially_removable = a.2;
            }

            // M Step
            // Find the best (most likely) count assignment
            let eq_vec = &self.eqb.count_vec;

            let (sender, receiver) = channel();

            eq_vec.par_iter().for_each_with(sender, |s, eqc| {
                let tg = &eqc.0;
                let v = &eqc.1;
                let csize = v.get_weights().len();
                let mut tmp_read_prob = vec![0.0; csize];
                let mut denom = 0.0;
                for read_mapping_cntr in 0..csize { //iterate through the list of strains in this equivalence class
                    let tgt = tg.get_tgts()[read_mapping_cntr];
                    if strain_valid[&tgt] {              // if the strain is valid, update the temp probability (score * current abundance * coverage 
                        let val = v.get_weights()[read_mapping_cntr] * strain_cnt[tgt] * self.strain_coverage[&tgt];
                        tmp_read_prob[read_mapping_cntr] = val;
                        denom += tmp_read_prob[read_mapping_cntr];
                    }
                }
                for read_mapping_cntr in 0..csize {
                    let tgt = tg.get_tgts()[read_mapping_cntr];
                    if strain_valid[&tgt] {
                        s.send((tgt, v.get_count() as f32 * (tmp_read_prob[read_mapping_cntr] / denom))).ok();
                    }
                }
            });

            //collect the information
            let res: Vec<_> = receiver.iter().collect();
            res.iter().for_each(|x| {
                let val = {
                    if x.1.is_nan() { 0.0 } else { x.1 }
                };
                new_strain_cnt[x.0] += val // adding as it is another
            });

            // E Step
            // normalize strain probabilities using the denum : p(s) = (count(s)/total_read_cnt)
            converged = true;
            let mut max_diff = 0.0;

            for i in 0..new_strain_cnt.len() {
                let adiff = (new_strain_cnt[i] - strain_cnt[i]).abs();
                if adiff > eps {
                    converged = false;
                }
                if adiff > max_diff {
                    max_diff = adiff;
                }
                strain_cnt[i] = new_strain_cnt[i];
                new_strain_cnt[i] = 0.0;
            }
            cntr += 1;
        }
        // input results into the strain_abundance variable. 
        let mut output_map = HashMap::new();
        let mut final_read_cnt = 0.0;
        let mut num_of_valids = 0;

        let valid_strains_cnt:Vec<f32> = strain_cnt.par_iter().enumerate().filter_map(|(index, val)| {
            if index == 0 { None }
            else if strain_valid[&index] {Some(*val)}
            else{ None }}).collect();

        let sum:f32 = valid_strains_cnt.par_iter().sum();

        for i in self.strain_abundance.keys() {
            final_read_cnt += strain_cnt[*i];
            if strain_valid[&i] {
                output_map.insert(self.ref_id_to_tax_id[i], strain_cnt[*i] / sum);
                num_of_valids += 1;
                // update taxa_abundance
                self.taxa_abundance.insert(self.ref_id_to_tax_id[i], strain_cnt[*i]);
            } else {
                output_map.insert(self.ref_id_to_tax_id[i], 0.0);
                self.taxa_abundance.insert(self.ref_id_to_tax_id[i], 0.0);
            }
        }
        println!("Final Reference-level read cnt: {}, # of valid refs: {}", final_read_cnt, num_of_valids);
        
        std::mem::swap(&mut self.strain_abundance, &mut output_map);
    }

    // outputs file with the references and their estimated abundancies
    pub(crate) fn serialize_simple(&mut self, output_filename: String) {
        println!("Write results into the file: {}", &output_filename);
        println!("# of strains: {}", self.strain_abundance.len());

        let mut output = File::create(output_filename).unwrap();
        for i in 1..self.strain_abundance.len() + 1 {
            let mut data = i.to_string();
            data.push_str("\t");
            data.push_str(&self.references.get(&i).unwrap().ref_name);
            data.push_str("\t"); 
            data.push_str(&self.strain_abundance[&i].to_string());
            data.push_str("\n");
            output.write(data.as_bytes()).ok();
        }
        println!("File has been written");
    }

    pub(crate) fn get_queries(&self) -> HashMap<usize, Query> {
        self.queries.clone()
    }

    pub(crate) fn get_references(&self) -> HashMap<usize, Reference> {
        self.references.clone()
    }

    pub(crate) fn get_strain_abundance(&self) -> HashMap<usize, f32> {
        self.strain_abundance.clone()
    }

    pub(crate) fn get_unmapping_reads(&self) -> usize {
        self.unmapping_reads
    }

    // main function for Cedar
    pub(crate) fn run_parallel(&mut self, mapper_output_name: String, 
        max_iter: usize,
        eps: f32,
        min_cnt: f32,
        segment_size: usize,
        flat_abundance: bool,
        batch_size: usize,
        method: String) {
            self.load_mapping_info_parallel(mapper_output_name, flat_abundance, segment_size, batch_size, method);
            self.parallel_em(max_iter, eps, min_cnt);
    }
}
