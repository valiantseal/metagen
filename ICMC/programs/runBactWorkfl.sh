
bash -i ./programs/main/download.sh

Rscript --vanilla ./programs/tests/checkDownloads.R

Rscript --vanilla ./programs/main/makeBacteriaLists.R

bash -i ./programs/main/forBactopia.sh

Rscript --vanilla ./programs/tests/bactopDataTest.R

bash -i ./programs/main/bactopia.sh

bash -i ./programs/main/pangenome.sh

bash -i ./programs/main/mlst.sh

Rscript --vanilla ./programs/main/iqtree.R

bash -i ./programs/main/plasmid.sh

# summaries
sh ./programs/main/sumAMR.sh

bash -i ./programs/main/sumBactopia.sh

Rscript --vanilla ./programs/main/sumMlst.R

sh ./programs/main/transferS3Custom.sh

# transfer custom output
# clean not needed directories
# transfer all output

bash -i ./programs/development/kleborate.sh