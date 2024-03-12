cd input

parallel gzip -d ::: *.fastq.gz
ls *_R1_*.fastq | parallel --verbose "sed -i 's# 1:N:#/1 N:#' {}"
ls *_R2_*.fastq | parallel "sed -i 's# 2:N:#/2 N:#' {}"
parallel gzip ::: *.fastq

# Create a list of unique sample prefixes by cutting off the filenames at the '_R1' or '_R2' part
for i in *R1_001.fastq.gz; do
    echo "${i%%_R1_001.fastq.gz}" >> ../newdir.list
done

cd ..

# Sort and remove duplicates
sort -u newdir.list > sorted_newdir.list
mv sorted_newdir.list newdir.list

# Create directories for each sample, and copy the corresponding 'R1' and 'R2' files into them
# Run clumpify and trimmomatic
for sample in $(cat newdir.list); do
    mkdir -p ./process/"$sample"/
    cp ./input/"${sample}_R1_001.fastq.gz" ./process/"$sample"/
    cp ./input/"${sample}_R2_001.fastq.gz" ./process/"$sample"/
    echo "$sample" > ./process/"$sample"/sample.name

    R1=./process/"$sample"/"${sample}_R1_001.fastq.gz"
    R2=./process/"$sample"/"${sample}_R2_001.fastq.gz"
    R1_clumped=./process/"$sample"/"${sample}_R1_001_clumped.fastq.gz"
    R2_clumped=./process/"$sample"/"${sample}_R2_001_clumped.fastq.gz"

    ~/bbmap/clumpify.sh in1=$R1 in2=$R2 out1=$R1_clumped out2=$R2_clumped dedupe

    R1_paired=./process/"$sample"/"${sample}_R1_001_paired.fastq.gz"
    R1_unpaired=./process/"$sample"/"${sample}_R1_001_unpaired.fastq.gz"
    R2_paired=./process/"$sample"/"${sample}_R2_001_paired.fastq.gz"
    R2_unpaired=./process/"$sample"/"${sample}_R2_001_unpaired.fastq.gz"

    java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $R1_clumped $R2_clumped $R1_paired $R1_unpaired $R2_paired $R2_unpaired ILLUMINACLIP:/home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    # Remove intermediate and unpaired files
    rm $R1 $R2 $R1_clumped $R2_clumped $R1_unpaired $R2_unpaired

    # Rename paired files back to original naming convention
    mv $R1_paired ./process/"$sample"/"${sample}_R1_001.fastq.gz"
    mv $R2_paired ./process/"$sample"/"${sample}_R2_001.fastq.gz"
done
