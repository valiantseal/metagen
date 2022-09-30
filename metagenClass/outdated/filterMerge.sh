conda activate ivar

sample=$(cat sample.name)
r1="$sample"_R1_001.fastq.gz
r2="$sample"_R2_001.fastq.gz

fastp -i "$r1" -I "$r2" \
--detect_adapter_for_pe \
--adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa \
--merge --merged_out merged_reads.fq.gz \
--out1 R1_not_merged.fq.gz \
--out2 R2_not_merged.fq.gz \
--thread 8