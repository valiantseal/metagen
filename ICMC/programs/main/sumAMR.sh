#AMR
mkdir antimic_res

# copy individual microbial resistances for each sample
for j in $(cat bacteria.list)
do
mkdir -p ./bactopia_gtdbtk/"$j"/antimic_res/combined
for i in ./bactopia_gtdbtk/"$j"/output/*/antimicrobial-resistance/*-report.txt; do cp "$i" ./bactopia_gtdbtk/"$j"/antimic_res/; done
done

echo ''
echo 'copied individual microbial resistances for each sample'
echo ''

# combine per bacteria antimicrobila resistances
for j in $(cat bacteria.list)
do
cd ./bactopia_gtdbtk/"$j"/antimic_res
Rscript --vanilla ../../../programs/main/combineAMR.R
cd ../../../
done

echo ''
echo 'combined per bacteria antimicrobila resistances'
echo ''

# copy combined for each species
cp ./bactopia_gtdbtk/*/antimic_res/combined/*.tsv ./antimic_res

# combine all bacteria that were used in a run together

cd antimic_res
mkdir combined
Rscript --vanilla ../programs/main/combAmrAllBact.R

echo ''
echo 'combined all bacteria that were used in a run together'
echo ''

cd ../
mkdir -p ./custom_output/antimicrobial_resistance
cp antimic_res/combined/* ./custom_output/antimicrobial_resistance