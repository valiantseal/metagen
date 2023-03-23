cd references
bowtie2-build MN908947.3.fna MN908947
cd ../

cd input
fastp -i EHC-C19-2884X_S16_L001_R1_001.fastq.gz -I EHC-C19-2884X_S16_L001_R2_001.fastq.gz \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w 8 -j fastp.json \
      -l 25

bowtie2 -p 8 -x ../references/MN908947 -1 filtered_1.fastq -2 filtered_2.fastq -S output.sam