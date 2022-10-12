conda activate seqtk

seqtk seq -a out.extendedFrags.fastq > merged_reads.fa
seqtk seq -a out.notCombined_1.fastq  > R1_not_merged.fa
seqtk seq -a out.notCombined_2.fastq > R2_not_merged.fa

echo merged_reads >> merged.list
echo R1_not_merged >> merged.list
echo R2_not_merged >> merged.list

for i in $(cat merged.list)
do
mkdir "$i"
mv "$i".fa ./"$i"/
echo "$i" > ./"$i"/merge.type
done
