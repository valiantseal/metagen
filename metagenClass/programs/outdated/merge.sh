sample=$(cat sample.name)
r1="$sample"_R1_001.fastq.gz
r2="$sample"_R2_001.fastq.gz

conda activate flash2
flash2 "$r1" "$r2" -t 8
