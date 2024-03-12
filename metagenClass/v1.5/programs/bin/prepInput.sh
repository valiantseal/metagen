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

# Make sure process_sample.sh is executable
chmod +x ~/v1.5/programs/bin/process_sample.sh

# Process each sample
cat newdir.list | parallel ~/v1.5/programs/bin/process_sample.sh {}
