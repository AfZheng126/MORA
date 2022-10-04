use std::collections::{HashMap, HashSet};

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
    position: i32,
    pub(crate) paired: bool
}

impl Mapping {
    fn new(reference_id: i32, score: usize, position: i32, paired: bool) -> Mapping {
        Mapping { reference_id, score, position, paired}
    }

    pub(crate) fn get_score(&self) -> f32 { self.score as f32}

    pub(crate) fn get_position(&self) -> i32 { self.position }

    pub(crate) fn get_reference_id(&self) -> usize {
        if self.reference_id != -1 { self.reference_id as usize} 
        else { usize::MAX }
    }
}

#[derive(Clone, PartialEq)]
pub(crate) struct Query {
    pub(crate) query_name: String, 
    cnt: usize,
    len: u32, 
    pub(crate) mappings: HashSet<Mapping>,        
    is_paired: bool, 
    total_score: f32,
}

impl Query {
    fn new(query_name: String, cnt: usize, len: u32, mappings: HashSet<Mapping>, is_paired: bool) -> Query {
        let mut total_score = 0.0;
        for map in &mappings {
            total_score += map.get_score();
        }
        Query { query_name, cnt, len, mappings, is_paired, total_score}
    }

    pub(crate) fn get_cnt(&self) -> usize { self.cnt }

    pub(crate) fn get_total_score(&self) -> f32 { self.total_score }

    fn add_mapping(&mut self, mapping: Mapping) {
        if mapping.reference_id != -1 {
            self.mappings.insert(mapping);
            self.cnt += 1;
            self.total_score += mapping.get_score();
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
use std::fs::File;
use std::io::{BufReader, BufRead};

use super::taxa::TaxaNode;

pub(crate) fn read_by_bash(file_name: String, batch_size: usize, method: String) -> (HashMap<usize, Reference>, HashMap<String, Query>){
    let f = File::open(file_name).unwrap();
    let stream = BufReader::with_capacity(batch_size, f);

    let (references, reference_to_id, stream, first_read) = analyze_header(stream);
    println!("references are done: {}", references.len());

    let queries = analyze_mappings_from_sam(stream, &reference_to_id, method, first_read);
    println!("queries are done: {}", queries.len());
    (references, queries)
}

fn analyze_header<R: BufRead>(mut stream: R) -> (HashMap<usize, Reference>, HashMap<String, usize>, R, String) {    
    let mut buffer = String::new();
    // get the header part of the SAM file
    let mut first_read = String::new();
    loop {
        let mut line = String::new();
        if stream.read_line(&mut line).unwrap() == 0 {
            break;
        }
        if line.starts_with("@") == false {
            first_read = line;
            // this should be the first read
            break;
        } else if line.starts_with("@SQ") {
            buffer = buffer + &line;
        }
    }
    let lines: Vec<&str> = buffer.split('\n').collect();

    // construct the references
    let references = (0..(lines.len()-1)).into_par_iter().fold(|| HashMap::new(), |mut acc, index| { //the -1 is because the last element of lines is ""
        let mut ref_name = String::new();
        let mut ref_len = 0;
        let chunks:Vec<_> = lines[index].split('\t').collect();
        // look at first tag
        let tag = chunks[0].as_bytes();
        if tag.len() !=3 || tag[0] != b'@' {
            panic!("Invalid header tag: {}\n", chunks[0]);
        }
        if tag == b"@SQ" { // if it is a reference, add it to the HashMap
            // look at the rest of the tags
            for i in 1..chunks.len() {
                let tag = chunks[i].as_bytes();
                if tag.len() < 3 || tag[2] != b':' {
                    panic!("Invalid tag: {} in line {}", chunks[i], i);
                }
                let tag_value = String::from_utf8(tag[3..].to_vec()).expect("Tag value not in UTF-8");
                if &tag[0..=1] == b"SN"{
                    ref_name = tag_value;
                } else if &tag[0..=1] == b"LN" {
                    ref_len = tag_value.parse().unwrap();
                }
            }
            acc.entry(index + 1).or_insert(Reference::new(ref_len, ref_name.to_string()));
        }
        acc
    }).reduce(|| HashMap::new(), |m1, m2| {
        m2.iter().fold(m1, |mut acc, (key, value)| {
            acc.entry(key.clone()).or_insert(value.clone());
            acc
        })
    });

    //duplicate the references to get a HashMap with keys and values switched. This is useful for files with a lot of references
    let reference_to_id = references.clone().into_par_iter().fold(|| HashMap::new(), |mut acc, (index, reference)| {
        acc.entry(reference.ref_name).or_insert(index);
        acc
    }).reduce(|| HashMap::new(), |m1, m2| {
        m2.iter().fold(m1, |mut acc, (key, value)| {
            acc.entry(key.clone()).or_insert(value.clone());
            acc
        })
    });
    (references, reference_to_id, stream, first_read)
}

fn line_to_mapping(line: Vec<&str>, references: &HashMap<String, usize>, method: &String) -> (String, bool, usize, Mapping) {
    let record_name = line[0].to_string();
    let flag_value:u32 = line[1].parse().unwrap();
    let is_paired = flag_value % 2 != 0;
    let record_length = line[9].len();
    if line[2] == "*" {
        let mapping = Mapping::new(-1, 0, 0, false);
        return (record_name, is_paired, record_length, mapping)
    }
    let reference_id = *references.get(&line[2].to_string()).unwrap() as i32;
    let position = line[3].parse().unwrap();
    let mut score:i32 = 0;

    //iterate through the tags to find the score
    for i in 11..line.len() {
        let tag = line[i].as_bytes();
        if tag.len() < 3 || tag[2] != b':' {
            panic!("Invalid tag: {} in line {}", line[i], i);
        }
        let mut tag_value = String::from_utf8(tag[5..].to_vec()).expect("Tag value not in UTF-8");
        tag_value = tag_value.replace("\n", "");
        if tag[0..=1] == vec![b'A', b'S'] {
            score = tag_value.parse().unwrap(); // it should be a positive integer, sometimes it was min(i32) for C++, but I considered that as 0.
            if method == "bowtie2" {
                score += 160;
            }
            if score <= 0 { // I assume that score = 0 means that it doesn't map
                score = 1
            }
        }
    }
    let mapping = Mapping::new(reference_id as i32, score as usize, position, is_paired);
    (record_name, is_paired, record_length, mapping) 
}

fn analyze_mappings_from_sam<R: BufRead>(mut stream: R, references: &HashMap<String, usize>, method: String, first_line: String) -> HashMap<String, Query> {
    let mut buffer = String::new();
    let mut queries:HashMap<String, Query> = HashMap::new();
    let mut left_over: String = "".to_string();
    let mut has_left_over = false;

    //work on first line left over from analyzing the header
    let line: Vec<_> = first_line.split('\t').collect();
    let (record_name, is_paired, record_length, mapping) = line_to_mapping(line, references, &method);
    let mut current_query = Query::new(record_name.to_string(), 0, record_length as u32, HashSet::new(), is_paired);
    current_query.add_mapping(mapping);
    // insert it as the current queries HashMap is empty
    queries.insert(record_name.to_string(), current_query);

    loop {
        buffer.clear();
        let buffer = std::str::from_utf8(stream.fill_buf().unwrap()).unwrap().to_string();
        let length = buffer.len();
        stream.consume(length);
        if length == 0 && !has_left_over{ break; }
        
        let buffer_string;
        if &left_over == "" {
            buffer_string = buffer;
        } else {
            buffer_string = [left_over.clone(), buffer].concat();
        }

        let mut lines: Vec<_> = buffer_string.split('\n').collect();

        //check if the last line of lines is a complete line
        left_over = lines.pop().unwrap().to_string();
        has_left_over = !(&left_over == "");
        if lines.len() == 0 { continue; }

        let line: Vec<_> = lines[0].split('\t').collect();
        let (record_name, is_paired, record_length, mapping) = line_to_mapping(line, references, &method);
        let mut current_query = Query::new(record_name.to_string(), 0, record_length as u32, HashSet::new(), is_paired);
        current_query.add_mapping(mapping);
        //check the first line to see if that query is already in our HashMap
        let mut current_name = record_name.clone();
        if queries.contains_key(&current_name) {
            for map in queries[&current_name].mappings.clone(){
                current_query.add_mapping(map);
            }
        }
        //go through the other lines in order
        for i in 1..lines.len() {
            let line: Vec<_> = lines[i].split('\t').collect();
            let (record_name, is_paired, record_length, mapping) = line_to_mapping(line, references, &method);
            if record_name == current_name {
                current_query.add_mapping(mapping);
            } else {
                queries.insert(current_name, current_query);
                current_name = record_name.clone();
                current_query = Query::new(record_name.to_string(), 0, record_length as u32, HashSet::new(), is_paired);
                current_query.add_mapping(mapping);
            }
        }
        queries.insert(current_name, current_query);
    }
    queries
}


pub(crate) fn read_rank(ref_2_tax_id_filename: String, batch_size: usize) -> HashMap<String, usize> {    
    let mut ref_2_tax_id = HashMap::new();
    let f = File::open(ref_2_tax_id_filename).unwrap();
    let mut stream = BufReader::with_capacity(batch_size, f);
    
    let mut left_over: String = "".to_string();
    let mut has_left_over = false;

    loop {
        let buffer = std::str::from_utf8(stream.fill_buf().unwrap()).unwrap().to_string();
        let length = buffer.len();
        stream.consume(length);
        if length == 0 && !has_left_over{ break; }
        
        let buffer_string;
        if &left_over == "" {
            buffer_string = buffer;
        } else {
            buffer_string = [left_over.clone(), buffer].concat();
        }

        let mut lines: Vec<_> = buffer_string.split('\n').collect();
        left_over = lines.pop().unwrap().to_string();
        has_left_over = !(&left_over == "");
        if lines.len() == 0 { continue; }

        let temp = (0..(lines.len())).into_par_iter().fold(|| HashMap::new(), |mut acc, index| {
            let chunks:Vec<_> = lines[index].split('\t').collect(); // chunks should be [name, id]
            let id: usize = chunks[1].parse().unwrap();
            acc.insert(chunks[0].to_string(), id);
            acc
        }).reduce(|| HashMap::new(), |m1, m2| {
            m2.iter().fold(m1, |mut acc, (key, value)| {
                // all the strings in the rank 2 taxId should be unique
                acc.entry(key.clone()).or_insert(value.clone());
                acc
            })
        });
        ref_2_tax_id.par_extend(temp);
    }
    ref_2_tax_id
}

pub(crate) fn read_taxonomy(taxonomy_tree_filename: String, batch_size: usize) -> HashMap<usize, TaxaNode>{
    let mut taxa_node_map = HashMap::new();
    
    let f = File::open(taxonomy_tree_filename).unwrap();
    let mut stream = BufReader::with_capacity(batch_size, f);
    
    let mut left_over: String = "".to_string();
    let mut has_left_over = false;

    loop {
        let buffer = std::str::from_utf8(stream.fill_buf().unwrap()).unwrap().to_string();
        let length = buffer.len();
        stream.consume(length);
        if length == 0 && !has_left_over{ break; }
        
        let buffer_string;
        if &left_over == "" {
            buffer_string = buffer;
        } else {
            buffer_string = [left_over.clone(), buffer].concat();
        }

        let mut lines: Vec<_> = buffer_string.split('\n').collect();
        left_over = lines.pop().unwrap().to_string();
        has_left_over = !(&left_over == "");
        if lines.len() == 0 { continue; }
        
        let temp = lines.into_par_iter().fold(|| HashMap::new(), |mut acc, line| {
            let chunks:Vec<_> = line.split('|').collect();
            if chunks.len() <= 2{
                acc
            } else {
                let id:usize = chunks[0].trim().parse().unwrap();
                let pid:usize = chunks[1].trim().parse().unwrap();
                let rank = chunks[2..].concat().trim().to_string();
                
                acc.insert(id, TaxaNode::new(id, pid, TaxaNode::str_to_rank(&rank)));

                acc
            }
        }).reduce(|| HashMap::new(), |m1, m2| {
            m2.iter().fold(m1, |mut acc, (key, value)| {
                // all the strings in the rank 2 taxId should be unique
                acc.entry(key.clone()).or_insert(value.clone());
                acc
            })
        });
        taxa_node_map.par_extend(temp);
    }
    taxa_node_map
}
