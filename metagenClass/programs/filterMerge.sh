conda activate ivar

for i in $(cat newdir.list);
do
cd process/"$i"

sample=$(cat sample.name)
r1="$sample"_R1_001.fastq.gz
r2="$sample"_R2_001.fastq.gz

fastp -i "$r1" -I "$r2" \
--detect_adapter_for_pe \
--adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa \
--merge --merged_out merged_reads.fastq \
--include_unmerged
--thread 8

cd ../../

done