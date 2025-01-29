# Bastilevirinae introns
Scripts used to analyse the data in the "Self-splicing introns in genes of <i>Bastillevirinae</i> bacteriophages" paper.
The analysis consists of four main steps:
1. Finding potential group I and group II introns in the genomes using Infernal matches to custom intron CM databases and resolving potential intron-exon boundaries based on alignment to known protein families.
2. Confirmation of introns splicing using cDNA reads and finding potential splicing sites.
3. Analysis of introns evolutionary relationships using cmalign and IQ-tree.
4. Preliminary prediction of secondary structure using LinearTurboFold and R2DT.


# REQUIREMENTS
Linux system (tested on Ubuntu 22.04 LTS)
conda or compatible package manager (dependencies are listed in the environment yaml and are installed automatically when the environment is created)


# INSTALLATION
Clone the repository and enter its main directory:
```
git clone https://github.com/zwmuam/Bastilevirinae_introns
cd Bastilevirinae_introns
```


Create environment based on environment.yml:
```
conda env create -f phage_introns_env.yml
```
OR
```
mamba env create -f phage_introns_env.yml
```


Before using you have to activate the environment:
```
conda activate phage_introns
```


Run the download_databases.sh script:
```
./download_databases.sh
```
By default the script downloads the original databases used in the paper to default locations used by the downstream tools.






# USAGE
Run the "find_introns" script to locate potential group I and group II matches:
```
./find_introns.py -f genomes.fasta -o genomes_intron_pred
```
This command assumes that default databases are installed but custom versions can be used using "-h" and "-c" arguments.
Optionally you may add a taxonomy xlsx table with "-w" to perform taxonomic distribution analysis.


If you have spliced cDNA reads similar to these described in the paper pass them to "confirm_introns" script to map transcripts to genome and identifies potential splicing sites:
```
./confirm_introns.py -ca adapter_promoter -cb barcodes_r -fq 100k_sample.fastq -rd rnaseq_references -o confirm_introns_test -s '__' -mr 1000
```
The script uses minimap2 to map reads to the genome and samtools to extract splicing sites. The "-s" argument is used to specify the separator between the barcode and the read name in the fastq file. The "-mr" argument is used to specify the maximum number of reads to process.


Introns for the alignment can be extracted using "get_intron_sequences":
```
./get_intron_sequences.py -f genomes.fasta -o introns.fna
```
By default the introns are extracted with 15-nt exon flanks compared to the original annotation to include any splicing-related elements like IGS etc but this can be changed using "-x" argument.


To re-annotate introns (useful if you added introns from external databases like GISSD) and generate R2DT-compatible annotations use "reannotate_introns":
```
./reannotate_introns.py -f GIISSD_and_Bastille_introns.fasta -o Annotated_GIISSD_and_Bastille_introns
```




To align all analysed introns we used InfeRNAl cmalign with "--matchonly" option (to exclude variable insertion sequences and focus on conserved RNA structures:
```
./cmalign -o Unique_vs_RF00028.sto --dnaout --mxsize 5140 --cpu 20 --matchonly --outformat Stockholm RF00028.cm GIISSD_and_Bastille.unique_099.fasta
```


The resulting stockholm alignment was converted to aligned fasta, columns with more than 50% gaps were masked in Geneious and run "prune_alignment" script to remove poorly aligned sequences:
```
./sto2fasta.py -s Unique_vs_RF00028.sto
# Mask columns in Genious
./prune_alignment.py -f Unique_vs_RF00028.fasta -l 50
```


Statistics for all unaligned sequences and dereplicate sequences as well as pruned and unpruned alignments were calculated using "alignment_statistics":
```
./alignment_statistics.py -u -f GIISSD_and_Bastille.fasta -o redundant_dataset_stats.xlsx
./alignment_statistics.py -u -f GIISSD_and_Bastille.unique_099.fasta -o nr_dataset_stats.xlsx
./alignment_statistics.py -f Unique_vs_RF00028.fasta -r 251 -o unpruned_aln_stats.xlsx
./alignment_statistics.py -f Unique_vs_RF00028.pruned.fasta -r 251 -o pruned_aln_stats.xlsx
```


The phylogenetic tree was constructed using IQ-tree:
```
iqtree -nt AUTO -s Unique_vs_RF00028.pruned.fasta -alrt 1000 -bb 1000
```
This command uses the automatic model selection and performs 1000 ultrafast bootstrap replicates and 1000 SH-aLRT replicates.


The secondary structure of the introns was predicted using LinearTurboFold:
```
./LinearTurboFold.py -f Unique_vs_RF00028.pruned.fasta -o Unique_vs_RF00028.pruned
```


LinearTurboFold alignements/structures, tsv annotations and json colour code (provided by reannotate_introns.py) were used to generate SVG visualisations using R2DT:
```
r2dt_from_gff.py -i All_introns_annotation.gff -o structure_annotation
# In the Docker container:
r2dt.py templatefree ${intron_structure}.fasta ${intron_structure}
python3 enrich_json.py --input-json ${intron_structure}.json --input-data ${intron_structure}.annot --output ${intron_structure}.enriched.json
json2svg.py -p colour_dict.json -i ${intron_structure}.enriched.json -o ${intron_structure}.enriched.svg
```
The visualisation is not included in the repository because it requires a fully configured R2DT pipeline.
We have run it within the R2DT Docker container based on the image provided at https://hub.docker.com/r/rnacentral/r2dt.
