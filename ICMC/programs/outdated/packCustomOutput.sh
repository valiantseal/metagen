mkdir custom_output
# pangenome
for i in $(cat bacteria.list)
do
mkdir custom_output/"$i"
cp ./bactopia_gtdbtk/"$i"/panOut/bactopia-tools/pangenome/pangenome/core-genome.distance.tsv ./custom_output/"$i"/
cp ./bactopia_gtdbtk/"$i"/panOut/bactopia-tools/pangenome/pangenome/iqtree/core-genome.contree ./custom_output/"$i"/
done

#AMR
mkdir antimic_res

# copy individual microbial resistances for each sample
for j in $(cat bacteria.list)
do
mkdir -p ./bactopia_gtdbtk/"$j"/antimic_res/combined
for i in ./bactopia_gtdbtk/"$j"/output/*/antimicrobial-resistance/*-report.txt; do cp "$i" ./bactopia_gtdbtk/"$j"/antimic_res/; done
done

# combine per bacteria antimicrobila resistances
for j in $(cat bacteria.list)
do
cd ./bactopia_gtdbtk/"$j"/antimic_res
Rscript --vanilla ../../../programs/combineRes.R
cd ../../../
done

# copy combined for each species
cp ./bactopia_gtdbtk/*/antimic_res/combined/*.tsv ./antimic_res

# combine all bacteria that were used in a run together

cd antimic_res
mkdir combined
Rscript --vanilla ../programs/combResAllBact.R
cd ../../
mkdir ./custom_output/antimicrobial_resistance
cp antimic_res/combined/* ./custom_output/antimicrobial_resistance

# combine assembled genomes from bactopia
mkdir ./custom_output/assembled_genomes
cp ./bactopia_gtdbtk/*/output/*/assembly/*.fna.gz ./custom_output/assembled_genomes/

# transfer
mkdir run_info
pwd > ./run_info/run.path
cat ./run_info/run.path |cut -f5 -d"/" > run_info/client.name
cat ./run_info/run.path |cut -f6 -d"/" > run_info/run.name

client=$(cat run_info/client.name)
run=$(cat run_info/run.name)

aws s3 cp --recursive ./custom_output/ s3://transfer-files-emory/ICMC/"$client"/"$run"/custom_output/
