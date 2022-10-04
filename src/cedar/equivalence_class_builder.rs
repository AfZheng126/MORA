use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use rayon::prelude::*;

//Hash function
fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

// maybe make TargetGroup and TGValue a single

#[derive(Eq, Hash, Clone)]
pub(crate) struct TargetGroup{
    /*
    tgts: A vector of potential targets in the references
    hash: the hash value?
    total_mass: ????
    valid: if the TargetGroup is a valid target group 
    */
    tgts: Vec<usize>,
    hash: u64,
    // total_mass: f32, 
    valid: bool,
}

impl TargetGroup{
    //The different constructors for TargetGroup
    fn new() -> TargetGroup {
        TargetGroup{
            tgts: Vec::new(), hash: 0, valid: true
        }
    }
    pub(crate) fn new_with_tgt(tgts: Vec<usize>) -> TargetGroup {
        let hash_value = calculate_hash(&tgts);
        TargetGroup{
            tgts, hash: hash_value, valid: true
        }
    }

    fn new_with_copy(other: &TargetGroup) -> TargetGroup {
        let mut tgts = Vec::new();
        for it in &other.tgts {
            tgts.push(it.to_owned());
        }
        TargetGroup {
            tgts, hash: other.hash, valid: other.valid}
    }
    
    // getters for attributes 
    pub fn get_tgts(&self) -> &Vec<usize> {
        &self.tgts
    }
    pub fn get_hash(&self) -> u64 {
        self.hash
    }
}

impl PartialEq for TargetGroup{
    fn eq(&self, other: &Self) -> bool {
        self.tgts == other.tgts  
    }
}
#[derive(Clone)]
pub(crate) struct TGValue {
    /*
    weights: vector over the references in the equivalence class, with value
        sum over reads of [score of read to this reference] / self.references[mapping.get_reference_id()].ref_len as f32
    combined_weights: 
    count: number of references in this equivalence class
    */
    weights: Vec<f32>,
    combined_weights: Vec<f32>,
    count: usize,
}

impl TGValue {
    //The constructors for TGValue
    fn new() -> TGValue {
        TGValue { weights: Vec::new(), combined_weights: Vec::new(), count: 0}
    }
    fn new_with_tgvalue_count(weights: Vec<f32>, count: usize) -> TGValue {
        TGValue { weights, combined_weights: Vec::new(), count}
    }
    fn new_from(other: &TGValue) -> TGValue {
        let mut weights = Vec::new();
        for it in &other.weights {
            weights.push(it.to_owned());
        }
        let mut combined_weights = Vec::new();
        for it in &other.combined_weights {
            combined_weights.push(it.to_owned());
        }

        TGValue { weights, combined_weights, count: other.count }
    }

    //normalizes the weights in self.weights
    fn normalize_aux(&mut self) {
        let mut sum_of_aux:f32 = 0.0;
        for i in 0..self.weights.len() {
            sum_of_aux += self.weights[i];
        }
        let norm: f32 = 1.0/sum_of_aux;
        for i in 0..self.weights.len() {
            self.weights[i] *= norm;
        }
    }

    pub fn get_count(&self) -> usize {
        self.count
    }

    pub fn get_weights(&self) -> &Vec<f32> {
        &self.weights
    }
}

#[derive(Clone)]
pub(crate) struct EquivalenceClassBuilder {
    //active_: whether the EquivalenceClassBuilder is active or not
    //count_map: a vector containing the TargetGroups and their TGValues
    //count_vec: a vector containing TargetGroups and their TGValues
    active_: bool,
    pub(crate) count_map: HashMap<TargetGroup, TGValue>,
    pub(crate) count_vec: Vec<(TargetGroup, TGValue)>, 
}

impl EquivalenceClassBuilder {

    pub(crate) fn new() -> EquivalenceClassBuilder {
        EquivalenceClassBuilder { active_: false, count_map: HashMap::new(), count_vec: Vec::new() }
    }

    //ends the EquivalenceClassBuilder
    pub(crate) fn finish(&mut self) {
        self.active_ = false;
        let mut total_count = 0;
        self.count_vec.reserve(self.count_map.len());
        for kv in &self.count_map {
            let mut new_tg_value = TGValue::new_from(kv.1);
            new_tg_value.normalize_aux();
            total_count += kv.1.count;
            self.count_vec.push((TargetGroup::new_with_copy(kv.0), new_tg_value));
        }
        println!("Counted {} total reads in the equivalence classes", total_count);
    }

    // Adds a group to the count_map with the weights if the target group is not in count_map,
    // otherwise, it adds the weights to the weights of the existing entry. 
    pub(crate) fn add_group(&mut self, g: TargetGroup, weights: Vec<f32>) {
        let index = self.count_vec.par_iter().position_any(|x| x.0 == g);
        if index == None {
            let v = TGValue::new_with_tgvalue_count(weights, 1);
            self.count_map.insert(g, v);
        } else {
            self.count_vec[index.unwrap()].1.count += 1;
            for i in 0..self.count_vec[index.unwrap()].1.weights.len() {
                self.count_vec[index.unwrap()].1.weights[i] += weights[i];
            }
        }
    }

    // merges another equivalenceClassBuilder's data into the current equivlanceClassBuilder's data
    // this does not delete the other EquivalenceClassBuilder
    fn merge_unfinished_eqb(&mut self, eqb: &EquivalenceClassBuilder) {
        for kv in eqb.count_map.keys() {
            let value = eqb.count_map.get(kv).unwrap();
            let new_key = TargetGroup::new_with_copy(kv);
            let new_value = TGValue::new_from(value);

            let it = self.count_map.entry(new_key).or_insert(new_value);
            for i in 0..it.weights.len() {
                it.weights[i] += value.weights[i];
            }
            it.count += value.count;
        }
    }
}