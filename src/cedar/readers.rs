use std::collections::{HashMap, HashSet};
use rust_htslib::bam::{Reader, Read,  Header, record::Aux};

use std::str;

#[derive(Clone, PartialEq, Eq, Hash, Copy)]
pub(crate) struct Mapping {
    /*
    stores the information of a mapping of the read to a reference (the read is known by which Query struct the Mapping is stored in)

    reference_id: the ID of the reference that this read is mapped to (-1 means no reference)
    score: the score of the mapping obtained from the SAM file
    position: the start position of the reference where the read is mapped to. 
    
    Currently not considering directionality
    */
    reference_id: i32,
    score: usize,
    position: i64,
    pub(crate) paired: bool
}

impl Mapping {
    fn new(reference_id: i32, score: usize, position: i64, paired: bool) -> Mapping {
        Mapping { reference_id, score, position, paired}
    }

    pub(crate) fn get_score(&self) -> f32 { self.score as f32}

    pub(crate) fn get_position(&self) -> i64 { self.position }

    pub(crate) fn get_reference_id(&self) -> usize {
        if self.reference_id != -1 { self.reference_id as usize} 
        else { usize::MAX }
    }
}

#[derive(Clone, PartialEq)]
pub(crate) struct Query {
    /*
    Query = read, otherwise there are too many things that start with an r. 

    query_id: the id to get the String of the query name
    cnt: the number of reference it maps to
    len: the length of the query
    total_score: sum of all the mapping scores
    unmapped: if cnt > 0
    */
    pub(crate) query_id: usize, 
    cnt: usize,
    len: u32, 
    pub(crate) mappings: HashSet<Mapping>,        
    is_paired: bool, 
    total_score: f32,
    unmapped: bool,
}

impl Query {
    fn new(query_name: usize, cnt: usize, len: u32, mappings: HashSet<Mapping>, is_paired: bool) -> Query {
        let mut total_score = 0.0;
        for map in &mappings {
            total_score += map.get_score();
        }
        let unmapped = {cnt != 0};
        Query { query_id: query_name, cnt, len, mappings, is_paired, total_score, unmapped}
    }

    pub(crate) fn get_cnt(&self) -> usize { self.cnt }

    pub(crate) fn get_total_score(&self) -> f32 { self.total_score }

    fn add_mapping(&mut self, mapping: Mapping) {
        if mapping.reference_id != -1 {
            self.mappings.insert(mapping);
            self.cnt += 1;
            self.total_score += mapping.get_score();
            self.unmapped = true;
        }
    }

    // returns the id and score of the best and second best mappings
    pub(crate) fn get_best_mappings(&self) -> (usize, f32, usize, f32) {
        let (mut highest_score, mut best_ref_id, mut second_highest_score, mut second_best_ref_id) = (0.0, 0, 0.0, 0);

        for mapping in &self.mappings {
            if mapping.get_score() > highest_score {
                best_ref_id = mapping.get_reference_id();
                highest_score = mapping.get_score();
            } else if mapping.get_score() > second_highest_score {
                second_best_ref_id = mapping.get_reference_id();
                second_highest_score = mapping.get_score();
            }
        }

        (best_ref_id, highest_score, second_best_ref_id, second_highest_score)
    }

    pub(crate) fn sort_mappings(&self) -> Vec<&Mapping> {
        let mut ordered_vec: Vec<&Mapping> = self.mappings.par_iter().collect();
        ordered_vec.par_sort_by(|a, b| b.get_score().partial_cmp(&a.get_score()).unwrap());
        ordered_vec
    }
}

#[derive(Clone)]
pub(crate) struct Reference {
    pub(crate) ref_len: usize,
    pub(crate) ref_name: String,
}

impl Reference {
    fn new(ref_len: usize, ref_name: String)-> Reference {
        Reference {ref_len, ref_name }
    }
}

use rayon::prelude::*;

/*
Inputs:
file_name: directory for SAM/BAM file
method: how to analyze mapping scores

Output:
(references, ref_id_2_names, queries, query_id_2_name)
*/
pub(crate) fn read_initial_alignments(file_name: String, method: String) -> (HashMap<usize, Reference>, HashMap<usize, Query>, HashMap<usize, String>){
    let f = Reader::from_path(&file_name).unwrap();
    let header = Header::from_template(f.header());
    let references = analyze_header(header);

    println!("references are done: {}", references.len());

    let (queries, query_id_2_name) = analyze_alignments(f, method);
    println!("reads are done: {}", queries.len());
    (references, queries, query_id_2_name)
}

fn analyze_header(header: Header) -> HashMap<usize, Reference> {
    let mut references = HashMap::new();
    for (key, records) in header.to_hashmap() {
        if key != "SQ".to_string() {
            continue;
        }
        references = records.par_iter().enumerate().fold(|| HashMap::new(), |mut acc, (ref_id, ref_record)| {
            let ref_name = &ref_record["SN"];
            let ref_len = ref_record["LN"].parse().unwrap();
            acc.entry(ref_id).or_insert(Reference::new(ref_len, ref_name.to_string()));
            acc
        }).reduce(|| HashMap::new(), |m1, m2| {
            m2.iter().fold(m1, |mut acc, (key, value)| {
                acc.entry(key.clone()).or_insert(value.clone());
                acc
            })
        });
    }
    references
}

/*
Inputs:
f: the thing that reads the SAM file
method: to determine how to analyze the mapping scores (mainly for bowtie2)
Output: 
(queries, names for queries)
*/
fn analyze_alignments(mut f: Reader, method: String) -> (HashMap<usize, Query>, HashMap<usize, String>) {
    let mut queries = HashMap::new();
    let mut query_id_2_names = HashMap::new();
    let mut query_id = 0;
    let mut current_name = "".to_string();
    let mut current_query = Query::new(0, 0, 0, HashSet::new(), true);
    
    for r in f.records() {
        let record = r.unwrap(); 
        let reference_id = record.tid();
        let query_name = str::from_utf8(record.qname()).unwrap();
        let mut score:i32 = 0;
        let position = record.pos();
        let paired = record.is_paired();
        let query_len = record.seq_len() as u32;
        match record.aux(b"AS") {
            Ok(value) => {
                if let Aux::U8(v) = value { score = v as i32; }
                if let Aux::I8(v) = value { score = v as i32; }
            } Err(_e) => {
                score = 0;
            }
        }
        if method == "bowtie2" { // This is where you add different ways to interperet AS:i scores for different intial aligners
            score += 160;
        }
        if score <= 0 { // I assume that score = 0 means that it doesn't map
            score = 1
        }
        let mapping = Mapping::new(reference_id, score as usize, position, paired);
        if query_name == current_name {
            current_query.add_mapping(mapping);
        } else {
            if query_id != 0 {
                queries.insert(query_id, current_query);
                query_id_2_names.insert(query_id, current_name.to_string());
            }
            query_id += 1;
            current_query = Query::new(query_id, 0, query_len, HashSet::new(), paired);
            current_name = query_name.to_string();
            current_query.add_mapping(mapping);
        }
    }    

    (queries, query_id_2_names)
}
