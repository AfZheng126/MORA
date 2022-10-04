use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

// SetCover Algorithm (Greedy Approach)
/*
sets: The potential covers
weigths: The weight of each element in the covers
unique_elements: The number of unique elements
*/
pub(crate) fn greedy_set_cover(sets: Vec<Vec<u64>>, weights: Vec<Vec<f32>>, uniqu_elements: usize) -> Vec<usize> {
    let mut tmp = HashSet::new();                     // what has currently been covered
    let mut left_over:HashSet<u64> = HashSet::new();                // what has not been covered 
    let mut final_cover = Vec::new();           // the output
    let mut covers = HashSet::new();                // ennumeration of the sets
    
    for i in 1..(uniqu_elements + 1) {
        left_over.insert(i as u64);
    }
    for i in 0..sets.len() {
        covers.insert(i);
    }
    
    while left_over.len() != 0  && covers.len() != 0 {
        // find the scores of the covers and add the lowest 5% of them to the final cover if it has at least 1 new element
        // this is done to decrease the amount to iterations needed to get a final cover
        // create the bins for the scores of the covers
        let mut score_bins: HashMap<usize, Vec<usize>> = covers.par_iter().fold(|| HashMap::new(), |mut acc, index| {
            let score = cover_score(&tmp, &sets[*index], &weights[*index]);
            let bin_num = (score * 2.0).floor() as usize;
            let entry = acc.entry(bin_num).or_insert(Vec::new());
            entry.push(*index);
            acc
        }).reduce(|| HashMap::new(), |m1, m2| {
            m2.iter().fold(m1, |mut acc, (key, value)| {
                let entry = acc.entry(key.clone()).or_insert(Vec::new());
                entry.extend(value);
                acc
            })
        });
        
        let mut cnt = 0;    // counter for how much elements has been added
        let mut bin = 0;    // bin score

        let mut loop_len = left_over.len() / 20;
        if loop_len < 1 {
            loop_len = 1;
        }
        while left_over.len() > 0 && cnt < loop_len {
            let mut bin_finished = true;
            if score_bins.contains_key(&bin) && score_bins[&bin].len() > 0 {
                while score_bins[&bin].len() > 0 {
                    // check conditions to continue
                    if cnt > loop_len || left_over.len() == 0 {
                        bin_finished = false;   // premature exit of the bin
                        break;
                    }
                    let mut can_push = false;
                    for k in &sets[score_bins[&bin][0]] {
                        if left_over.remove(k) {
                            tmp.insert(*k);
                            can_push = true;
                        }
                    }
                    if can_push {
                        final_cover.push(score_bins[&bin][0]);
                    }
                    let entry = score_bins.entry(bin).or_insert(Vec::new());
                    entry.remove(0);
                    cnt += 1;
                }
            }
            if bin_finished { bin += 1; }
        }

        // update covers
        covers.clear();
        for (_, v) in score_bins {
            for k in v {
                covers.insert(k);
            }
        }
    }
    final_cover
}

// Helper function for calculating the score of a potential cover. The score is calculated as 
// sum of the weights of each element / # of elements in s that is not in vs. 
fn cover_score(vs: &HashSet<u64>, s: &Vec<u64>, weights: &Vec<f32>) -> f32{
    let mut v = 0.0;
    for i in 0.. s.len() {
        v += weights[i];
    }
    let mut count = 0.0; // keeps track of the number of elements in s that is not in vs. 
    for element in s {
        if !vs.contains(&element){
            count += 1.0;
        }
    }
    v/count
}
