# create database directory
mkdir -p ./databases
# Download databases from the paper (PHROGs v4 and custom intron-related InfeRNAl models).
wget -O ./databases/databases.zip "https://figshare.com/ndownloader/files/47854990"
unzip ./databases/databases.zip
hmmpress ./databases/Phrogs4_HMMer3.hmm
cmpress ./databases/Merged.1.GISSD_IRFAM.cm


# If you want to replace PHROGs with a newer version, you need to recalculate HMMer models (see commented code below)
# Exact URLs and paths may vary from version to version
: '
wget -O ./databases/phrog_annot_v4.tsv "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv" # update this to download newer version
wget -O ./databases/MSA_phrogs.tar.gz "https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz" # update this to download newer version
tar -xvzf databases/MSA_phrogs.tar.gz -C databases/
for a in ./databases/MSA_Phrogs_M50_FASTA/*.fma; do hmmbuild "${a}.hmm" "$a"; done
for h in ./databases/MSA_Phrogs_M50_FASTA/*.hmm; do cat "$h" >> ./databases/meged.hmm; rm "$h"; done
hmmpress databases/MSA_phrogs.hmm
'
