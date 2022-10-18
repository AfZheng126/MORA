extern crate clap;
use clap::{Arg, App};
use std::{time::Instant, collections::HashMap, fs::File, io::{BufReader, BufRead}};

mod cedar;
use cedar::Cedar;

mod assignment;
use crate::assignment::{assign_mappings, write_output, write_output_with_taxonomy};

fn main() {
    let commands = App::new("LRA").version("1.0").author("126andrew.zheng@gmail.com")
                        .about("Long-Read Re-Alignment")
                        .arg(Arg::with_name("SAM File")
                            .short('s')
                            .long("sam")
                            .help("path for sam file")
                            .required(true)
                            .takes_value(true)
                            .display_order(1))
                        .arg(Arg::with_name("Min Count")
                            .short('c')
                            .long("minCnt")
                            .help("minimum count for a reference to be valid")
                            .takes_value(true)
                            .default_value("0.1")
                            .display_order(3))
                        .arg(Arg::with_name("Abund output")
                            .short('a')
                            .long("abund_out")
                            .help("path for abundance output")
                            .takes_value(true)
                            .display_order(4))
                        .arg(Arg::with_name("Batch Size")
                            .short('b')
                            .long("batch_size")
                            .help("batch_size for reading SAM")
                            .default_value("100000000")
                            .takes_value(true)
                            .display_order(3))
                        .arg(Arg::with_name("Max EM iterations")
                            .long("max_em")
                            .help("maximum allowed iterations of EM")
                            .takes_value(true)
                            .default_value("300")
                            .display_order(3))
                        .arg(Arg::with_name("segment size")
                            .long("segment_size")
                            .help("size to split referenes into")
                            .takes_value(true)
                            .default_value("100")
                            .display_order(3))
                        .arg(Arg::with_name("Min score diff")
                            .long("min_score_diff")
                            .help("minimum difference between mapping scores divided by best mapping score for second step")
                            .takes_value(true)
                            .default_value("0.5")
                            .display_order(3))
                        .arg(Arg::with_name("Max abund diff")
                            .long("max_abund_diff")
                            .help("maximum difference between abundance levels")
                            .takes_value(true)
                            .default_value("0.001")
                            .display_order(3))
                        .arg(Arg::with_name("Output")
                            .short('o')
                            .long("output")
                            .help("path for final output of assignments")
                            .required(true)
                            .takes_value(true)
                            .display_order(2))
                        .arg(Arg::with_name("taxonomy")
                            .long("tax")
                            .help("write output with taxonomy details with provided tax directory")
                            .takes_value(true)
                            .display_order(4))
                        .arg(Arg::with_name("Method")
                            .long("method")
                            .help("mapping method: (pufferfish, bowtie2, minimap2)")
                            .takes_value(true)
                            .display_order(2)
                            .default_value("pufferfish"))
                        .arg(Arg::with_name("Threads")
                            .short('t')
                            .long("threads")
                            .help("number of threads for rayon to use")
                            .takes_value(true)
                            .display_order(4)
                            .default_value("3"))
                        .get_matches();

    // collect values from user inputs
    let sam_file = commands.value_of("SAM File").unwrap().to_string();
    let min_cnt: f32 = commands.value_of("Min Count").unwrap().parse().unwrap();
    let max_iter: usize = commands.value_of("Max EM iterations").unwrap().parse().unwrap();
    let batch_size: usize = commands.value_of("Batch Size").unwrap().parse().unwrap();
    let segment_size: usize = commands.value_of("segment size").unwrap().parse().unwrap();
    let method: String = commands.value_of("Method").unwrap().to_string();
    let threads:usize = commands.value_of("Threads").unwrap().parse().unwrap();
    
    // setup number of threads
    rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();

    let mut cedar = Cedar::new(); 
    cedar.run_parallel(sam_file, max_iter, 0.001, min_cnt, segment_size, batch_size, method);

    if commands.is_present("Abund output") {
        cedar.serialize_simple(commands.value_of("Abund output").unwrap().to_string());
    }

    let abund_diff: f32 = commands.value_of("Max abund diff").unwrap().parse().unwrap();
    let score_diff: f32 = commands.value_of("Min score diff").unwrap().parse().unwrap();
    let output_filename = commands.value_of("Output").unwrap();

    let output = assign_mappings(cedar, abund_diff, score_diff);

    println!("\nWriting results to {}", &output_filename);
    if commands.is_present("taxonomy") {
        let tax_dir = commands.value_of("taxonomy").unwrap().to_string();
        println!("tax directory = {}", &tax_dir);
        write_output_with_taxonomy(output_filename.to_string(), output, 
                                        tax_dir.to_string() + "/accessionsTaxIDs.tab",
                                        tax_dir.to_string() + "/nodes.dmp",
                                        tax_dir.to_string() + "/names.dmp");
    } else {
        write_output(output_filename.to_string(), output);
    }
}
