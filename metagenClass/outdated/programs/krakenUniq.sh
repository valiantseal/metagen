conda activate krakenuniq

for j in $(cat merged.list)
do 
cd "$j"

sample=$(cat merge.type)

krakenuniq --threads 8 --db /home/ubuntu/extraVol/krakenUniq/viral \
--report-file krakUniq_sample.report --output krakUniq_sample.kraken \
--classified-out krakUniq_classified.reads "$sample".fa

cd ../
done
