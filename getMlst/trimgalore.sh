r1=$(cat ./sample.name)_R1_001.fastq.gz
r2=$(cat ./sample.name)_R2_001.fastq.gz


/home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20  --paired "$r1" "$r2" \
--basename "trimmed" --cores 2

