# MORA

Mora is an aligner that re-aligns query reads to a unique reference.

Main steps of Mora: 
1. Find the potential mappings of the query reads and output them as a SAM file using an initial aligner.
2. Calculate the expected abundance levels of the references based on the SAM file.
3. Assign each query that had at least one valid mapping to a reference based on their mapping scores and the expected abundance levels.
4. Output the results into a txt file. 

# Requirements
[Rust](https://www.rust-lang.org/tools/install) and Cargo need to be installed and added to PATH.

# Installation
```
git clone https://github.com/AfZheng126/MORA.git
cd MORA
bash install.sh
cargo build --release
```

# Running MORA
After everything in the config files (see below) is updated according to your directories, run 
```
snakemake --snakefile MORA --cores 24 --resources mem_mb=140000
```

# Running MORA as a Rust Program
If you already have a SAM file that has mappings scores stored in the AS:i: optional field, you can directly run the Rust program and skip the indexing and mapping steps. To do this and get outputs without taxonomic information, run the following commands in the Mora directory.
```
cargo run --release -- -s samfile -o output
```
For more options and customization, run 
```
cargo run -- -h
```
If you are runing from another directory and the specific binary is wanted, run 
```
target/release/mora -h
```

# Config File
The parameters of the config.yaml file used for the snakemake pipline are listed below: 

| Parameter | Description |
| ---- | --- |
| BINARIES | Binary folder directory (default: binaries) - do not edit |
| REFERENCES | Directory to reference fasta file |
| SAMPLES_DIR | Directory to folder containing query fasta files |
| RESULTS | Directory to write the results |
| FILES_EXT | Query files extension, i.e. .fq, .fq.gz etc |
| MAPPING_MODE | Algorithm for the initial mapping - (pufferfish, bowtie2, minimap2)|
| STRATEGY | "PE" for paired-end samples or "SE" for single-end samples |
| TYPE | RNA or DNA host-specific samples - right now only supports DNA |
| BATCH_SIZE | Size of buffer for reading sam file |
| MIN_CNT | Minimum number of counts for a reference to be considered valid |
| MIN_SCORE_DIFFERENCE | Minimum score difference for a query to be assgined second |
| MAX_ABUNDANCE_DIFFERENCE | Maximum difference allowed between the initial abundance estimation and the abundances created from assignments |
| SEGMENT_SIZE | Size to split references into bins |
| ABUNDANCE_OUTPUT | Whether to output estimated abundance levels |
| TAXONOMY | Directory of taxonomic information to write results with taxonomic classes (NA to not include taxonomic information in the results) |
| MEM_MB | Amount of memory to be allocated to snakemake |
| TPS | Number of threads to be used per sample |

# Query Files
The program requires a list of query files. These can be .fasta, .fq, or even compressed files. If the query files are pair-end quries, their name must be of the form *_1.fq and *_2.fq, where the file extension can something else. The directory of these query files must be written into the config file. 

# Reference File
If a reference file is provided, its directory must also be written into the config file. If there is no reference file, you can download the fasta file representing the complete representative and reference bacterial genomes from NCBI RefSeq database by following the instructions from the [Microbial reference preparation](https://github.com/ivlachos/agamemnon/wiki/Use-case) from the Agamemnon Wiki. The index will be built when you run the program, so you don't have to manually do it. 

# Taxonomic Information
Normally, the output of the program is two columns telling you which reference each query came from. If taxonomic information about the assigned reference is wanted as well, extra files must be made. To do this, navigate to the scripts directory and run
```
bash taxonomy.sh reference.fa
```
where reference.fa is your reference files. After this is done, update TAXONOMY in the config file to Taxonomy. 


 # Use Case
Sample data and the results from the Mora paper can be found [here](https://github.com/AfZheng126/MORA-data). To run the data, simply update the configuration file with where you download the data and run it with the snakemake command. 
