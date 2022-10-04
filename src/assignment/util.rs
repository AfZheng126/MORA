fn assignments_2_lineage(assignments: HashMap<String, String>, lineage: HashMap<String, Vec<Lineage>>, accessions_2_tax: HashMap<String, String>) -> HashMap<String, (String, Vec<Lineage>, Vec<Lineage>)> {
    let mut assignments2lineage = HashMap::new();
    for (query, reference) in &assignments {
        let mut assigned_lineage = create_na_lineage();
        let mut real_lineage = create_na_lineage();
        let chunks: Vec<_> = query.split(".").collect();
        let query_name = chunks[0].to_string();
        if accessions_2_tax.contains_key(reference) {
            assigned_lineage = lineage[&accessions_2_tax[reference]].clone();
        }
        if accessions_2_tax.contains_key(&query_name) {
            real_lineage = lineage[&accessions_2_tax[&query_name]].clone();
        }
    
        assignments2lineage.insert(query.to_string(),(reference.to_string(), assigned_lineage, real_lineage));
    }
    assignments2lineage
}