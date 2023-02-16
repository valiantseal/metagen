# too many intermidiate files
r1=$(ls *R1.fastq) 

samtools view -f 16 -b -o negative-strand.bam virema1.bam
samtools view -F 16 -b -o positive-strand.bam virema1.bam

machine=$(head -n 1 "$r1" | sed 's/:/\t/g' - | awk '{print $1}' -)
    
samtools bam2fq negative-strand.bam | grep "$machine" - > samp_negative-mapped-reads_temp.txt
samtools bam2fq positive-strand.bam | grep "$machine" - > samp_positive-mapped-reads_temp.txt
Rscript --vanilla ~/extraVol/Copyback/nextflowTrial/bin/getNegUniqName.R

negRead=$(ls *_negative-mapped-reads_unique_temp.txt) # edit to work on exact file names
mv "$negRead" samp_negative-mapped-reads_unique_temp.txt

sort -n samp_positive-mapped-reads_temp.txt | uniq -d - >> samp_positive-mapped-reads_dedup_temp.txt
sort -n samp_positive-mapped-reads_temp.txt | uniq -u - >> samp_positive-mapped-reads_dedup_temp.txt
sed "s/@//g" samp_positive-mapped-reads_dedup_temp.txt > samp_positive-mapped-reads_dedup_temp1.txt
sed "s/@//g" samp_negative-mapped-reads_unique_temp.txt > samp_negative-mapped-reads_unique_temp1.txt

seqtk subseq virema1.fastq samp_negative-mapped-reads_unique_temp1.txt > samp_neg-mapped-reads_temp.fastq
seqtk subseq virema1.fastq samp_positive-mapped-reads_dedup_temp1.txt > samp_pos-mapped-reads_temp.fastq
seqtk seq -r samp_neg-mapped-reads_temp.fastq > samp_neg-mapped-rev_temp.fastq # reverse sequence
    
cat samp_neg-mapped-rev_temp.fastq samp_pos-mapped-reads_temp.fastq > final_virema.fastq