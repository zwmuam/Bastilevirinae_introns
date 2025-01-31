# Bastilevirinae introns
Scripts used to analyse the data in the "Self-splicing introns in genes of <i>Bastillevirinae</i> bacteriophages" paper.
The analysis consists of four main steps:
1. Finding potential group I and group II introns in the genomes using Infernal matches to custom intron CM databases and resolving potential intron-exon boundaries based on alignment to known protein families.
2. Confirmation of introns splicing using cDNA reads and finding potential splicing sites.
3. Analysis of introns evolutionary relationships using cmalign and IQ-tree.
4. Preliminary prediction of secondary structure using LinearTurboFold and R2DT.


# REQUIREMENTS
Linux system (tested on Ubuntu 22.04 LTS), Conda or compatible package manager (dependencies are listed in the phage_introns_env.yml and installed automatically when the environment is created)


# INSTALLATION
Clone the repository and enter the main directory:
```
git clone https://github.com/zwmuam/Bastilevirinae_introns
cd Bastilevirinae_introns
```


Create the environment:
```
conda env create -f phage_introns_env.yml
```
OR
```
mamba env create -f phage_introns_env.yml
```


Activate the environment:
```
conda activate phage_introns
```


Run the download_databases.sh script:
```
./download_databases.sh
```
By default the script downloads the database versions used in the paper to default location (./databases/).


# USAGE
Locate initial group I and group II matches:
```
./find_introns.py -f genomes.fasta -o genomes_intron_pred
```
This command assumes that default databases are installed. Custom HMMs and CMs can be used "-h" and "-c" arguments.
Optionally, xlsx table with taxonomy can be added as "-w" argument to perform distribution analysis.


If you have spliced cDNA reads (see the paper) pass them to "confirm_introns" script to map transcripts to genome and identify splicing sites:
```
./confirm_introns.py -ca custom_adapters -cb custom_barcodes -fq reads.fastq -rd reference_genomes_dir -o output_dir -s '__' -mr 1000
```
The script uses modified [Porechop](https://github.com/rrwick/Porechop) to demultiplex reads and remove adapters. Clean reads are mapped to reference genomes with minimap2 and analysed by spliced_bam2gff to find splicing sites. The "-s" argument is used to specify the separator between the barcode and the read name in the fastq file. The "-mr" specifies the number of reads in bam files used for downstream visualiation


To align all analysed introns to the reference model from RFAM:
```
./cmalign -o cmalign_output.sto --matchonly --outformat Stockholm reference.cm non_redundant_seqences.fasta
```
We used [InfeRNAl](http://eddylab.org/infernal) cmalign with "--matchonly" option and RF00028 model to exclude variable insertion sequences and focus on conserved RNA structures. In case of our machine we also added "--mxsize 5140" "--cpu 20" parameters to increase the maximum matrix size as well as number of threads and improve the performance. Some tools (e.g. tree construction programs) may also require DNA-like output triggered by "--dnaout" option.


Convert stockholm ouput of cmalign to aligned fasta:
```
./sto2fasta.py -s cmalign_output.sto
```


Columns with more than 50% gaps can be masked using one of commonly used sequence analysis software suites. To prune poorly aligned sequences run:
```
./prune_alignment.py -f cmalign_output.fasta -l 50
```


To re-annotate introns (useful for intron sequences imported from external databases like GISSD) and generate R2DT-compatible annotations use:
```
./reannotate_introns.py -f redundant_sequences.fasta -o Annotated_redundant_sequences
```


Generate statistics for all sequence datasets and alignments:
```
./alignment_statistics.py -u -f redundant_sequences.fasta -o redundant_dataset_stats.xlsx
./alignment_statistics.py -u -f non_redundant_seqences.fasta -o nr_dataset_stats.xlsx
./alignment_statistics.py -f cmalign_output.fasta -r 251 -o unpruned_aln_stats.xlsx
./alignment_statistics.py -f cmalign_output.pruned.fasta -r 251 -o pruned_aln_stats.xlsx
```


Construct a phylogenetic tree:
```
iqtree -nt AUTO -s Unique_vs_RF00028.pruned.fasta -alrt 1000 -bb 1000
```
This command uses [IQ-tree](http://www.iqtree.org) with automatic model selection, 1000 ultrafast bootstrap replicates and 1000 SH-aLRT replicates.


Predict the secondary structure of introns from different clades observed in the phylogenetic tree:
```
for clade in clade_*.fna; do ./linearturbofold -i ${clade} -o ${clade}_ltf --pf -v --bpp; done
```
[LinearTurboFold](https://github.com/LinearFold/LinearTurboFold) is run to simultaneously align and fold the sequences from each clade ("-v" prints alignment, folding and runtime information, --pf saves partition function and "--bpp" writes base pair probabilities).
This step is not included in the environment. To get the software, please visit the original repository by Li and colleagues.


Use structures, tsv annotations and json colour code (provided by reannotate_introns.py) to generate structure visualisations using R2DT:
```
# In the Docker container:
r2dt.py templatefree ${intron_structure}.fasta ${intron_structure}
python3 enrich_json.py --input-json ${intron_structure}.json --input-data ${intron_structure}.annot --output ${intron_structure}.enriched.json
json2svg.py -p colour_dict.json -i ${intron_structure}.enriched.json -o ${intron_structure}.enriched.svg
```
The visualisation is not included in the repository because it requires a fully configured [R2DT pipeline](https://docs.r2dt.bio/en/latest/about.html).
Run it within the R2DT Docker container based on the [image from Docker Hub](https://hub.docker.com/r/rnacentral/r2dt).
