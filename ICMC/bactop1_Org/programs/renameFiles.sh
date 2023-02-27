cd input
for i in *R1_001.fastq.gz
do
newName=$(echo "$i" | cut -f1,2,3 -d "_").fastq.gz
mv "$i" "$newName"
done 

cd ../

cd input
for i in *R2_001.fastq.gz
do
newName=$(echo "$i" | cut -f1,2,3 -d "_").fastq.gz
mv "$i" "$newName"
done 

cd ../

bactopia prepare ./input > fastqs.txt