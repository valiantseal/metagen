
cd ./input

for i in *_R1_*
do
newName=$(echo "$i" | cut -f1 -d"_")_R1.fastq.gz
mv "$i" "$newName"
done


for i in *_R2_*
do
newName=$(echo "$i" | cut -f1 -d"_")_R2.fastq.gz
mv "$i" "$newName"
done