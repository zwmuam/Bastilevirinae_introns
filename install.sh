
# create conda environment based on environment.yml
conda env create -f phage_introns_env.yml -p ./venv

# download databases
# download PHROGs database and (from RFAM and GSIID)
wget -O ./databases/MSA_phrogs.tar.gz "https://phrogs.lmge.uca.fr/downloads_from_website/MSA_phrogs.tar.gz"
tar -xvzf databases/MSA_phrogs.tar.gz -C databases/
wget -O ./databases/phrog_annot_v4.tsv "https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv"

# convert databases to binary format
hmmpress databases/MSA_phrogs.hmm
cmpress databases/MSA_phrogs.cmp

echo "Installation complete. Please activate the conda environment with 'conda activate ./venv' before running the pipeline."
