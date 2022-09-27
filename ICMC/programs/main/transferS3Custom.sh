mkdir custom_output
# pangenome
for i in $(cat bacteria.list)
do
mkdir -p custom_output/pangenome/"$i"
cp ./bactopia_gtdbtk/"$i"/panOut/bactopia-tools/pangenome/pangenome/core-genome.distance.tsv ./custom_output/pangenome/"$i"/
cp ./bactopia_gtdbtk/"$i"/panOut/bactopia-tools/pangenome/pangenome/iqtree/core-genome.contree ./custom_output/pangenome/"$i"/
done

#AMR
mkdir ./custom_output/antimicrobial_resistance
cp antimic_res/combined/* ./custom_output/antimicrobial_resistance

# combine assembled genomes from bactopia
mkdir ./custom_output/assembled_genomes
cp ./bactopia_gtdbtk/*/output/*/assembly/*.fna.gz ./custom_output/assembled_genomes/

# bactopia summary
cp ./bactopia_summary/combined/summary_combined.csv ./custom_output/quality_bactopia_summary_combined.csv

#MLST
cp ./mlst_summary/mlst_types.txt ./custom_output/mlst_types.csv

# transfer
mkdir run_info
pwd > ./run_info/run.path
cat ./run_info/run.path |cut -f5 -d"/" > run_info/client.name
cat ./run_info/run.path |cut -f6 -d"/" > run_info/run.name

client=$(cat run_info/client.name)
run=$(cat run_info/run.name)

aws s3 cp --recursive ./custom_output/ s3://transfer-files-emory/ICMC/"$client"/"$run"/custom_output/

echo s3://transfer-files-emory/ICMC/"$client"/"$run"/ > download.me
aws s3 cp download.me s3://transfer-files-emory/