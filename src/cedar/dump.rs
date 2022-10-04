/* 
        fn process_reads(&mut self) -> Stats {
        //let (mut total_read_cnt, mut seq_not_found, mut total_multi_mapping_reads, mut total_unmapped_reads, 
        //   mut total_reads_not_passing_cond, mut tid, mut conflicting) = (0, 0, 0, 0, 0, 0, 0);
        let (mut total_read_cnt, mut total_multi_mapped_reads, mut total_unmapped_reads) = (0, 0, 0);

        for (_, query) in self.queries.clone() {

            //println!("current query: {}, cnt: {}, ", &self.queries[i].r_name, self.queries[i].get_cnt());

            if query.get_cnt() != 0 {
                if query.get_cnt() > 1 {
                    total_multi_mapped_reads += 1;
                }
                let mut mapping_scores = Vec::new();
                let mut read_per_strain_prob_inst = Vec::new();     // Vector of (reference ID, weighted mapping score)
                for mapping in &query.mappings {
                    let ref_name = &self.references[&mapping.get_reference_id()].ref_name;
                    // ignore reference that we don't have a taxaID for
                    if self.ref_name_2_tax_id.contains_key(ref_name) {
                        let entry = self.cov.entry(self.ref_name_2_tax_id[ref_name]).or_insert(0);
                        *entry += mapping.get_score() as usize;
                        mapping_scores.push((mapping.get_score(), mapping.paired));
                        read_per_strain_prob_inst.push((mapping.get_reference_id(), mapping.get_score() / self.references[&mapping.get_reference_id()].ref_len as f32));
                    }
                }
                // update the eqb
                self.update_eqb(read_per_strain_prob_inst);
            } else {
                total_unmapped_reads += 1;
            }
            total_read_cnt += query.get_cnt();
        }
        Stats::new_with_stats(total_read_cnt, total_multi_mapped_reads, total_unmapped_reads)
    }
    
    // update the EquivalenceClassBuilder
    fn update_eqb(&mut self, read_per_strain_prob_inst: Vec<(usize, f32)>) {
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
        self.eqb.add_group(tgt, probs);
    }

    pub(crate) fn load_mapping_info(&mut self, mapper_output_filename: String, require_concordance: bool, flat_abundance: bool, only_unique: bool, only_perfect: bool, segment_size: usize, range_factorization_bins: u32, batch_size: usize) {
        println!("Cedar: Load Mapping File");
        println!("Mapping Ouput File: {}", mapper_output_filename);
        
        // load the information from the file
        let c = read_by_bash(mapper_output_filename, batch_size);     // right now I am only considering unpaired reads
        self.references = c.0;
        self.queries = c.1;

        println!("Data from SAM file has been read");

        //initially guess equal abundance levels
        let value = 1.0/ self.references.len() as f32;
        for key in self.references.keys() {
            self.strain_abundance.insert(*key, value);
        }        
        println!("Constucting coverage bins");
        // construct converage bins (splits the references into subsequences of a given length)
        
        for i in 0..self.references.len() {
            let ref_name = &self.references.get(&(i + 1)).unwrap().ref_name;
            let ref_length = self.references.get(&i).unwrap().ref_len;
            let bin_cnt = ref_length / segment_size + 1;
            let bins = vec![0; bin_cnt];
            self.strain_coverage_bins.insert(i, bins);
            let tid;  // the ID of the taxonomy of the reference. 
            if flat_abundance {
                tid = i;
            } else {
                tid = self.ref_name_2_tax_id.get(ref_name).unwrap().to_owned();
            }
            self.ref_id_to_tax_id.insert(i, tid);
            self.cov.insert(tid, 0);
        }


        // update the information in the Cedar struct
        let stats = self.process_reads();
        self.read_cnt = stats.get_total_read_cnt();
        self.update_bins(segment_size);
        self.calculate_coverage();

        // print the information obtained from the mapping_output_file 
        stats.print_stats();
    }    

    pub(crate) fn basic_em(&mut self, max_iter: usize, eps: f32, min_cnt: f32) {
        self.eqb.finish();

        let eq_vec = &self.eqb.count_vec;

        let max_seq_id = self.strain_abundance.len();

        let mut new_strain_cnt = vec![0.0; max_seq_id];
        let mut strain_valid:HashMap<usize, bool> = HashMap::new();
        let mut strain_potentially_removable:HashMap<usize, bool> = HashMap::new();

        for i in 1..(max_seq_id + 1) {
            strain_valid.insert(i, true);
            strain_potentially_removable.insert(i, false);
        }
        

        let mut strain_cnt = vec![0.0; max_seq_id + 1];

        for (key, value) in &self.strain_abundance {
            strain_cnt[*key] = *value;
        }
        //println!("the updated strain_cnt is {:?}", &strain_cnt);

        // logger stuff
        println!("max_seq_Id = {}", max_seq_id);
        println!("found: {} equivalence classes", eq_vec.len());

        let mut tot_count:usize = 0;
        for eqc in eq_vec {
            tot_count += eqc.1.get_count();
        }
        // logger stuff
        println!("Total starting count {}", tot_count);
        println!("Total mapped reads cnt {}", self.read_cnt); 

        let mut cntr:usize = 0;
        let mut converged = false;
        let thresholding_iter_step = 10;
        let mut can_help = true;
        //async::threadpool_schedular asyncScheduler(numThreads) -- Later when incorporating threads

        // let block_size:u32 = 1000;

        // not sure what this does. It seems to create a vector of tuples but not sure how and what the name is
        //std::vector<std::pair<decltype(eq_vec.begin()), decltype(eqvec.begin())>>
        //range(std::ceil(static_cast<double >(eq_vec.size())/ static_cast<double>(blockSize)));
        
        // let idx:u64 = 0;
        // for i in (0..eq_vec.len()).step_by(block_size) {

        // }

        while cntr < max_iter && converged == false {
            if cntr % thresholding_iter_step == 0 && can_help {
                let a = self.apply_set_cover(&strain_cnt, strain_valid, strain_potentially_removable, min_cnt, can_help);
                println!("applied set cover");
                can_help = a.0;
                strain_valid = a.1;
                strain_potentially_removable = a.2;
            }


            println!("M step");
            // M Step
            // Find the best (most likely) count assignment
            // There is some stuff with async which I will not include for now. 
            let eq_vec = &self.eqb.count_vec;

            for eqc in eq_vec { // for each equivalence class, we update all the strain abundances in it
                let tg = &eqc.0;
                let v = &eqc.1;
                let csize = v.get_weights().len();
                let mut tmp_read_prob = vec![0.0; csize];
                let mut denom = 0.0;
                for read_mapping_cntr in 0..csize { //iterate through the list of strains in this equivalence class
                    let tgt = tg.get_tgts()[read_mapping_cntr];
                    if strain_valid[&tgt] {              // if the strain is valid, update the temp probability (score * current abundance * coverage 
                        tmp_read_prob[read_mapping_cntr] = v.get_weights()[read_mapping_cntr] * strain_cnt[tgt] * self.strain_coverage[&tgt];
                        denom += tmp_read_prob[read_mapping_cntr];
                    }
                }
                for read_mapping_cntr in 0..csize {
                    let tgt = tg.get_tgts()[read_mapping_cntr];
                    if strain_valid[&tgt] { 
                        new_strain_cnt[tgt] += v.get_count() as f32 * (tmp_read_prob[read_mapping_cntr] / denom) // adding as it is another 
                    }
                }
            }

            println!("E Step");
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
        println!("iterator cnt: {}", cntr);

        // input results into the strain_abundance variable. 
        let mut output_map = HashMap::new();
        output_map.reserve(self.strain_abundance.len());
        let mut final_read_cnt = 0.0;
        let mut num_of_valids = 0;

        let mut sum = 0.0;
        for i in 0..strain_cnt.len() {
            sum += strain_cnt[i];
        }

        for i in 0..self.strain_abundance.len() {
            final_read_cnt += strain_cnt[i];
            if strain_valid[&i] {
                num_of_valids += 1;
                if let Some(ele) = output_map.get_mut(&i) {
                    *ele += strain_cnt[i] / sum;
                }
                // update taxa_abundance
                self.taxa_abundance.insert(self.ref_id_to_tax_id[&(i+1)], strain_cnt[i]);
            } else {
                self.taxa_abundance.insert(self.ref_id_to_tax_id[&(i+1)], 0.0);
            }
        }

        println!("Final Reference-level read cnt: {}, # of valid refs: {}", final_read_cnt, num_of_valids);

        std::mem::swap(&mut self.strain_abundance, &mut output_map);
    }



    pub(crate) fn run_single(&mut self, mapper_output_name: String, 
        require_concordance: bool,
        max_iter: usize,
        eps: f32,
        min_cnt: f32,
        output_name: String,
        only_unique: bool,
        only_perfect: bool,
        segment_size: usize,
        range_factorization_bins: u32,
        flat_abundance: bool,
        batch_size: usize) {
            self.load_mapping_info(mapper_output_name, require_concordance, flat_abundance, only_unique, only_perfect, segment_size, 
                range_factorization_bins, batch_size);
            self.basic_em(max_iter, eps, min_cnt);
            if !flat_abundance {
                //self.serialize(output_name, "STRAIN".to_string());
                self.serialize_simple(output_name);
            } 
    }

        // writes a file containing the following information
    // taxa_id    taxa_rank    taxa_count   
    fn serialize(&mut self, output_filename: String, prune_level_in: String) {
        println!("Write results into the file: {}", &output_filename);
        println!("# of strains: {}", self.strain_abundance.len());

        let mut valid_taxa:HashMap<usize, f32> = HashMap::new();

        for i in 0..self.strain_abundance.len() {
            if self.taxa_node_map.contains_key(&i)  {
                let mut walker = self.taxa_node_map.get(&i).unwrap();
                while !walker.is_root() && walker.rank != TaxaNode::str_to_rank(&prune_level_in) {
                    walker = self.taxa_node_map.get(&(walker.get_parent_id() as usize)).unwrap();
                    if walker.get_id() > 18000000000000 {   // Not sure what decides the big number 
                        std::process::exit(1);              // Not sure which exit code to use here, maybe will create a custum one. 
                    }
                }
                
                if !walker.is_root() {
                    if valid_taxa.contains_key(&walker.get_id()) {
                        *valid_taxa.get_mut(&walker.get_id()).unwrap() += self.strain_abundance[&i];
                    } else {
                        valid_taxa.insert(walker.get_id(), self.strain_abundance[&i]);
                    }
                }
            } else {
                panic!("Taxa not found");
            }
        }

        let mut final_read_cnt = 0.0;
        let mut output = File::create(output_filename).unwrap();
        for kv in valid_taxa {
            if kv.1 != 0.0 {
                let mut data = kv.0.to_string();
                data.push_str("\t");
                data.push_str(&TaxaNode::rank_to_string(&self.taxa_node_map.get(&kv.0).unwrap().rank));
                data.push_str("\t");
                data.push_str(&kv.1.to_string());
                data.push_str("\n");
                output.write(data.as_bytes()).ok();
            }
            final_read_cnt += kv.1;
        }
        println!("Final reported read count by Cedar: {}", final_read_cnt);
        // We could also add the cov numbers. 
    }

*/