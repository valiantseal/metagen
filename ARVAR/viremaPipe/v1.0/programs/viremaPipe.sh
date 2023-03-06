r1=$(ls *R1.fastq) 
r2=$(ls *R2.fastq); 

fastp -i "$r1" -I "$r2" \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w 4 -j fastp.json \
      -m --merged_out merged_prep_temp.fastq -A -l 25 \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa
      

perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_1.fastq filtered_1_fastp-tagged_temp.fastq 1 & perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_2.fastq filtered_2_fastp-tagged_temp.fastq 2;
sh /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/editMerged.sh
cat filtered_1_fastp-tagged_temp.fastq filtered_2_fastp-tagged_temp.fastq merged_prep_temp-tag.fastq >virema1.fastq
    

python /home/ubuntu/extraVol/Copyback/test_builds_2885y/programs/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna virema1.fastq \
    virema1.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p 4 -Overwrite >> ViReMaR1-results.txt
    

samtools view -buSh -F 4 virema1.sam | samtools sort --threads 4 - -o virema1.bam
samtools view -f 16 -b -o negative-strand.bam virema1.bam
samtools view -F 16 -b -o positive-strand.bam virema1.bam

machine=$(head -n 1 "$r1" | sed 's/:/\t/g' - | awk '{print $1}' -)
    
samtools bam2fq negative-strand.bam | grep "$machine" - > samp_negative-mapped-reads_temp.txt
samtools bam2fq positive-strand.bam | grep "$machine" - > samp_positive-mapped-reads_temp.txt
Rscript --vanilla ~/extraVol/Copyback/nextflowTrial/bin/getNegUniqName.R

negRead=$(ls *_negative-mapped-reads_unique_temp.txt)
mv "$negRead" samp_negative-mapped-reads_unique_temp.txt

sort -n samp_positive-mapped-reads_temp.txt | uniq -d - >> samp_positive-mapped-reads_dedup_temp.txt
sort -n samp_positive-mapped-reads_temp.txt | uniq -u - >> samp_positive-mapped-reads_dedup_temp.txt
sed "s/@//g" samp_positive-mapped-reads_dedup_temp.txt > samp_positive-mapped-reads_dedup_temp1.txt
sed "s/@//g" samp_negative-mapped-reads_unique_temp.txt > samp_negative-mapped-reads_unique_temp1.txt

seqtk subseq virema1.fastq samp_negative-mapped-reads_unique_temp1.txt > samp_neg-mapped-reads_temp.fastq
seqtk subseq virema1.fastq samp_positive-mapped-reads_dedup_temp1.txt > samp_pos-mapped-reads_temp.fastq
seqtk seq -r samp_neg-mapped-reads_temp.fastq > samp_neg-mapped-rev_temp.fastq
    
cat samp_neg-mapped-rev_temp.fastq samp_pos-mapped-reads_temp.fastq > final_virema.fastq

python /home/ubuntu/extraVol/Copyback/test_builds_2885y/programs/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna final_virema.fastq \
    final_virema.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p 4 -Overwrite >> final_virema-results.txt