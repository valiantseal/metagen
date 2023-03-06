r1=$(ls *R1.fastq) 
r2=$(ls *R2.fastq); 

fastp -i "$r1" -I "$r2" \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w 16 -j fastp.json \
      -m --merged_out merged_prep_temp.fastq -A -l 25 \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
      

perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_1.fastq filtered_1_fastp-tagged_temp.fastq 1 & perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_2.fastq filtered_2_fastp-tagged_temp.fastq 2;
sh /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/editMerged.sh
cat filtered_1_fastp-tagged_temp.fastq filtered_2_fastp-tagged_temp.fastq merged_prep_temp-tag.fastq >virema1.fastq
