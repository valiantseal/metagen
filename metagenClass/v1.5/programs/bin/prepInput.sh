cd input

# Create a list of unique sample prefixes by cutting off the filenames at the '_R1' or '_R2' part
for i in *R1_001.fastq.gz; do 
    echo "${i%%_R1_001.fastq.gz}" >> ../newdir.list
done

cd ..

# Sort and remove duplicates
sort -u newdir.list > sorted_newdir.list
mv sorted_newdir.list newdir.list

# Create directories for each sample, and copy the corresponding 'R1' and 'R2' files into them
for sample in $(cat newdir.list); do
    mkdir -p ./process/"$sample"/
    cp ./input/"$sample"_R1_001.fastq.gz ./process/"$sample"/
    cp ./input/"$sample"_R2_001.fastq.gz ./process/"$sample"/
    echo "$sample" > ./process/"$sample"/sample.name
done
