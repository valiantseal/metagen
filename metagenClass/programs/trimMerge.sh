sample=$(cat sample.name)
r1="$sample"_R1_001.fastq.gz
r2="$sample"_R2_001.fastq.gz


#/home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 30 --length 36 --paired "$i"_1.fastq.gz "$i"_2.fastq.gz
/home/ubuntu/TrimGalore-0.6.6/trim_galore --quality 20 --paired "$r1" "$r2" --basename "trimmed" --cores 4


conda activate flash2
flash2 trimmed_val_1.fq.gz trimmed_val_2.fq.gz -t 12


