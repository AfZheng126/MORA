# Mora

Mora is an read re-assigner that re-assigns query reads to a unique reference.

Main steps of Mora: 
1. Calculate the expected abundance levels of the references based on the input SAM file.
3. Assign each query that had at least one valid mapping to a reference based on their mapping scores and the expected abundance levels.
4. Output the results into a txt file. 

For more details, please consult the (preprint) [paper](https://www.biorxiv.org/content/10.1101/2022.10.18.512733v2).

# Requirements
Rust and Cargo need to be installed to and added to PATH. The instructions can be found on the [Rust website](https://www.rust-lang.org/tools/install) or in the [Rust cookbook](https://doc.rust-lang.org/cargo/getting-started/installation.html) if you are using the command line.

# Installation

1. [rust](https://www.rust-lang.org/tools/install) **version > 1.60.0** and associated tools such as cargo are required and assumed to be in PATH.
2. Standard tools such as gcc, make, cmake may be needed. 

To install and download the necessary rust packages to run MORA, run the following commands. 
```
git clone https://github.com/AfZheng126/MORA.git
cd MORA
cargo build --release
```

# Running Mora as a Rust Program
If you already have a SAM file that has mappings scores stored in the AS:i: optional field, you can directly run the Rust program and skip the indexing and mapping steps. To do this and get outputs without taxonomic information, run the following commands in the Mora directory.
```
cargo run --release -- -s samfile -o output
```
If you are runing from another directory and the specific binary is wanted, run 
```
target/release/mora -s samifile -o output
```
A sample sam file is provided in the samples directory. To use it, use the following command. 
```
target/release/mora -s sample/test.sam -o output.txt
```
For more options and customization, run 
```
target/release/mora -h
```

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
