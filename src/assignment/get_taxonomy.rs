use std::{collections::{HashMap, HashSet}, fs::File, io::{BufReader, BufRead, Write}};

#[derive(Clone)]
struct Lineage {
    rank: String,
    name: String,
}

impl Lineage {
    fn new(rank: String, name: String) -> Lineage {
        Lineage { rank, name }
    }
    fn get_rank(&self) -> String {
        self.rank.to_string()
    }
    fn get_name(&self) -> String {
        self.name.to_string()
    }
}

/* 
    Add NA to the missing taxonomy ranks
*/
fn fix_empty_ranks(lineage: HashMap<String, Vec<Lineage>>) -> HashMap<String, Vec<Lineage>> {
    let ranks = vec!["species".to_string(), "genus".to_string(), "family".to_string(), "order".to_string(), "class".to_string(), "phylum".to_string(), "superkingdom".to_string()];
    let mut new_lineage = HashMap::new();
    for (key, mut val) in lineage {
        let mut new_val = Vec::new();
        // push first lineage as it is the key
        new_val.push(val[0].clone());
        val.remove(0);
        // fill the rest in terms of rank
        for rank in 0..ranks.len(){
            let mut new_element = Lineage::new(ranks[rank].to_string(), "NA".to_string());
            for i in 0..val.len() {
                if val[i].get_rank() == ranks[rank]{
                    new_element = val[i].clone();
                    val.remove(i);
                    break;
                }
            }
            new_val.push(new_element);
        }
        new_lineage.insert(key, new_val);
    }
    new_lineage
}

/*
    assignments: <query name, reference name>
    at_file: accessions 2 taxIDs file
    OUTPUT: <TaxID,  vector of accession numbers that have the given TaxID>
*/
fn temp_accessions_2_tax_id(assignments: &HashMap<String, String>, at_file: String) -> (HashMap<String, Vec<String>>, HashMap<String, String>) {
    let mut tax_id_accesions = HashMap::new();
    let mut accessions_2_tax = HashMap::new();
    let mut accessions = HashSet::new();
    for (_, reference) in assignments {
        accessions.insert(reference.to_string());
    }
    let f = File::open(at_file).unwrap();
    let stream = BufReader::new(f);
    let lines: Vec<_> = stream.lines().collect();
    
    for line in lines {
        let line = line.unwrap();
        let chunks:Vec<_> = line.split('\t').collect();
        let entry = tax_id_accesions.entry(chunks[1].to_string()).or_insert(Vec::new());
        entry.push(chunks[0].to_string());
        accessions_2_tax.insert(chunks[0].to_string(), chunks[1].to_string());
    }
    let accessions_found:HashSet<String> = tax_id_accesions.keys().cloned().collect();
    let mut not_found = HashSet::new();
    for element in accessions {
        if !accessions_found.contains(&element) {
            not_found.insert(element);
        }
    }
    println!("number of accessions not found: {}", not_found.len());
    (tax_id_accesions, accessions_2_tax)
}

/*
    nodes_file: nodes file
    names_file: names file
    tax_id_accessions: <TaxID,  vector of accession numbers that have the given TaxID>
    OUTPUT: 
*/
fn build_taxonomy(nodes_file: String, names_file: String, tax_id_accessions: &HashMap<String, Vec<String>>) -> HashMap<String, Vec<Lineage>>{
    let (mut nodes, mut names, mut lineage) = (HashMap::new(), HashMap::new(), HashMap::new());
    let ranks = vec!["species".to_string(), "genus".to_string(), "family".to_string(), "order".to_string(), "class".to_string(), "phylum".to_string(), "superkingdom".to_string()];
    
    let f = File::open(nodes_file).unwrap();
    let stream = BufReader::new(f);
    let lines: Vec<_> = stream.lines().collect();

    // nodes: <TaxID, <Parent TaxID, rank of TaxID>>
    // example: <10, <1706371, genus>>
    for line in lines {
        let line = line.unwrap();
        let chunks:Vec<_> = line.split('|').collect();
        nodes.insert(chunks[0].trim().to_string(), vec![chunks[1].trim().to_string(), chunks[2].trim().to_string()]);
    }
    
    // names: <TaxID, scientific name>
    let f = File::open(names_file).unwrap();
    let stream = BufReader::new(f);
    let lines: Vec<_> = stream.lines().collect();

    for line in lines {
        let line = line.unwrap();
        let chunks:Vec<_> = line.split('|').collect();
        if chunks[3].trim() == "scientific name" {
            names.insert(chunks[0].trim().to_string(), chunks[1].trim().to_string());
        }
    }
    
    for key in tax_id_accessions.keys() {       //key is TaxID
        let mut temp_lineage = Vec::new();  // rank| name of parent| parent TaxID
        let mut helpful_list = Vec::new();
        let mut parent = key;
        while parent != "131567" && parent != "10239" {
            let val = Lineage::new(nodes[parent][1].to_string(), names[parent].to_string());
            temp_lineage.push(val);
            parent = &nodes[parent][0];
        }
        let mut rank_counter = 0;
        for rank in 0..temp_lineage.len() {
            if rank_counter == 0 {
                helpful_list.push(temp_lineage[rank].clone());
            } else {
                let parent_rank = temp_lineage[rank].get_rank();
                if ranks.contains(&parent_rank) { // checks if the rank of parent is valid
                    helpful_list.push(temp_lineage[rank].clone());
                }
            }
            rank_counter += 1;
        }
        lineage.insert(key.to_owned(), helpful_list); // lineage: <TaxID, list of TaxID and its ancestors in the form: " rank of key | name of key"
    }
    lineage = fix_empty_ranks(lineage);
    lineage
}

// find the lineages of the assignments
fn assignments_2_lineage(assignments: HashMap<String, String>, lineage: HashMap<String, Vec<Lineage>>, accessions_2_tax: HashMap<String, String>) -> HashMap<String, (String, Vec<Lineage>)> {
    let mut assignments2lineage = HashMap::new();
    for (query, reference) in &assignments {
        let mut assigned_lineage = create_na_lineage();
        if accessions_2_tax.contains_key(reference) {
            assigned_lineage = lineage[&accessions_2_tax[reference]].clone();
        }    
        assignments2lineage.insert(query.to_string(),(reference.to_string(), assigned_lineage));
    }
    assignments2lineage
}


fn create_na_lineage() -> Vec<Lineage> {
    let mut v = Vec::new();
    for _ in 0..8{
        v.push(Lineage::new("NA".to_string(), "NA".to_string()));
    }
    v
}


/*
Write the results out in the following format:
Query reads     Assigned Reference      Species     Genus       Family      Order       Class       Phylum      Superkingdom

Inputs: 
    out_dir: where to write the results
    assignment2lineage: the results
*/

fn write_results(out_dir: String, assignment2lineage: HashMap<String, (String, Vec<Lineage>)>) {
    let mut results = Vec::new();
    let mut output = File::create(out_dir).unwrap();
    let data = ["Query", "Reference", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Superkingdom"].join("\t");
    output.write(data.as_bytes()).ok();
    output.write("\n\n".as_bytes()).ok();
    for (query, (reference, assigned)) in &assignment2lineage {
        let mut result = Vec::new();
        result.push(query.to_string());
        result.push(reference.to_string());
        for i in 1..assigned.len() {
            result.push(assigned[i].get_name());
        }
        results.push(result);
    }
    for result in results {
        let i = result.join("\t");
        output.write(i.as_bytes()).ok();
        output.write("\n".as_bytes()).ok();
    }
}


pub(crate) fn tax_main(assignments: HashMap<String, String>, at_file: String, nodes_file: String, names_file: String, out_dir: String) {
    let (tax_id_accesions, accessions_2_tax) = temp_accessions_2_tax_id(&assignments, at_file); // map the accessions to the taxonomic IDs
    let lineage = build_taxonomy(nodes_file, names_file, &tax_id_accesions); // create lineages of the references
    let assignments2lineage = assignments_2_lineage(assignments, lineage, accessions_2_tax); //
    write_results(out_dir, assignments2lineage);
}
