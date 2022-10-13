use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

//Hash function
fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

#[derive(Eq, Hash, Clone)]
pub(crate) struct TargetGroup{
    /*
    tgts: A vector of potential targets in the references
    hash: the hash value?
    valid: if the TargetGroup is a valid target group 
    */
    tgts: Vec<usize>,
    hash: u64,
    valid: bool,
}

impl TargetGroup{
    //The different constructors for TargetGroup
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
    count: number of reads in this equivalence class
    */
    weights: Vec<f32>,
    combined_weights: Vec<f32>,
    count: usize,
}

impl TGValue {
    //The constructors for TGValue
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
    //count_map: a HashMap containing the TargetGroups and their TGValues
    active_: bool,
    pub(crate) count_map: HashMap<TargetGroup, TGValue>,
}

impl EquivalenceClassBuilder {

    pub(crate) fn new() -> EquivalenceClassBuilder {
        EquivalenceClassBuilder { active_: false, count_map: HashMap::new() }
    }

    //ends the EquivalenceClassBuilder
    pub(crate) fn finish(&mut self) {
        self.active_ = false;
        let mut total_count = 0;
        let mut new_map = HashMap::new();
        for (tg, val) in &self.count_map {
            let mut new_tg_value = TGValue::new_from(val);
            new_tg_value.normalize_aux();
            total_count += val.count;
            new_map.insert(TargetGroup::new_with_copy(tg), new_tg_value);
        }
        self.count_map = new_map;
        println!("Counted {} total reads in {} equivalence classes", total_count, self.count_map.len());
    }

    // Adds a group to the count_map with the weights if the target group is not in count_map,
    // otherwise, it adds the weights to the weights of the existing entry. 
    pub(crate) fn add_group(&mut self, g: TargetGroup, weights: Vec<f32>) {
        if self.count_map.contains_key(&g) {
            let mut tg_val = TGValue::new_from(self.count_map.get(&g).unwrap());
            tg_val.count += 1;
            for i in 0..tg_val.weights.len() {
                tg_val.weights[i] += weights[i];
            }
            self.count_map.insert(g, tg_val);
        } else {
            let v = TGValue::new_with_tgvalue_count(weights, 1);
            self.count_map.insert(g, v);
        }
    }
}