conda activate seqtk

seqtk seq -a merged_reads.fq.gz > merged_reads.fa
seqtk seq -a R1_not_merged.fq.gz > R1_not_merged.fa
seqtk seq -a R2_not_merged.fq.gz > R2_not_merged.fa

conda activate kraken2

mkdir kraken2
for i in *.fa;
do
echo "$i" >> ../../kraken2.log
kraken2 --use-names --threads 8 --db /home/ubuntu/kraken2/viral \
--report ./kraken2/"$i".report  "$i" > ./kraken2/"$i".kraken
done


