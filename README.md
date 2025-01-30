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
This command assumes that default databases are installed. Custom HMMs and CMs can be used using "-h" and "-c" arguments.
Optionally, xlsx table with taxonomy can be addes as "-w" argument to perform distribution analysis.


If you have spliced cDNA reads (see the paper) pass them to "confirm_introns" script to map transcripts to genome and identify splicing sites:
```
./confirm_introns.py -ca adapter_promoter -cb barcodes_r -fq 100k_sample.fastq -rd rnaseq_references -o confirm_introns_test -s '__' -mr 1000
```
The script uses mofified porchop to demultiplex reads and remove adapters. Clean reads are mapped to reference genomes with minimap2 and analysed by spliced_bam2gff to find splicing sites. The "-s" argument is used to specify the separator between the barcode and the read name in the fastq file. The "-mr" specifies the number of reads in bam files used for downstream visualiation


To align all analysed introns to the referenc group I model from RFAM:
```
./cmalign -o Unique_vs_RF00028.sto --dnaout --mxsize 5140 --cpu 20 --matchonly --outformat Stockholm RF00028.cm GIISSD_and_Bastille.unique_099.fasta
```
We used InfeRNAl cmalign with "--matchonly" option to exclude variable insertion sequences and focus on conserved RNA structures.


Convert stockholm ouput of cmalign to aligned fasta:
```
./sto2fasta.py -s Unique_vs_RF00028.sto
```


Columns with more than 50% gaps can be masked using one of commonly used sequence analysis software suites. To prune poorly aligned sequences run:
```
./prune_alignment.py -f Unique_vs_RF00028.fasta -l 50
```


To re-annotate introns (useful for intron sequences imported from external databases like GISSD) and generate R2DT-compatible annotations use:
```
./reannotate_introns.py -f GIISSD_and_Bastille_introns.fasta -o Annotated_GIISSD_and_Bastille_introns
```


Generate statistics for all sequence datasets and alignments:
```
./alignment_statistics.py -u -f GIISSD_and_Bastille.fasta -o redundant_dataset_stats.xlsx
./alignment_statistics.py -u -f GIISSD_and_Bastille.unique_099.fasta -o nr_dataset_stats.xlsx
./alignment_statistics.py -f Unique_vs_RF00028.fasta -r 251 -o unpruned_aln_stats.xlsx
./alignment_statistics.py -f Unique_vs_RF00028.pruned.fasta -r 251 -o pruned_aln_stats.xlsx
```


Construct a phylogenetic tree:
```
iqtree -nt AUTO -s Unique_vs_RF00028.pruned.fasta -alrt 1000 -bb 1000
```
This command uses IQ-tree with automatic model selection, 1000 ultrafast bootstrap replicates and 1000 SH-aLRT replicates.


Predict the secondary structure of introns:
```
linearturbofold -i clade_8.fna -o clade_8_ltf  -v --pf --bpp
```
LinearTurboFold is run to XXX clades sim align


Use structures, tsv annotations and json colour code (provided by reannotate_introns.py) to generate structure visualisations using R2DT:
```
# In the Docker container:
r2dt.py templatefree ${intron_structure}.fasta ${intron_structure}
python3 enrich_json.py --input-json ${intron_structure}.json --input-data ${intron_structure}.annot --output ${intron_structure}.enriched.json
json2svg.py -p colour_dict.json -i ${intron_structure}.enriched.json -o ${intron_structure}.enriched.svg
```
The visualisation is not included in the repository because it requires a fully configured R2DT pipeline.
Run it within the R2DT Docker container based on the image provided at https://hub.docker.com/r/rnacentral/r2dt.
