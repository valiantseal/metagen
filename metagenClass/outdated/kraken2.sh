conda activate kraken2

for j in $(cat merged.list)
do 
cd "$j"
sample=$(cat merge.type)

kraken2 --use-names --threads 8 --db /home/ubuntu/kraken2/viral \
--report sample.report  "$sample".fa > sample.kraken
cd ../
done


