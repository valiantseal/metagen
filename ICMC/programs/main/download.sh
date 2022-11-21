conda activate fastq-dl

cat samples.list | parallel -j 8 'fastq-dl {} SRA -o ./input'