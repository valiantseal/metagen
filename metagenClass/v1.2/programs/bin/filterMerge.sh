cd process

ls -d */ | parallel -j 4 'cd {} && fastp -i *R1_001.fastq.gz -I *R2_001.fastq.gz \
--detect_adapter_for_pe \
--adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa \
--merge --merged_out merged_reads.fastq \
--include_unmerged \
-l 100 \
--thread 8'