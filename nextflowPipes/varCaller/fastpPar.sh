cd gnuPar
ls -d */ | parallel -j 6 'cd {} && r1=$(ls *R1.fastq); r2=$(ls *R2.fastq); 
fastp -i "$r1" -I "$r2" \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w 8 -j fastp.json \
      -m --merged_out merged_prep_temp.fastq -A -l 25 \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa' 