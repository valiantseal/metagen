conda activate fastq-dl2

mkdir -p input

cat sra_samples.list | parallel -j 7 'fastq-dl --accession {} --provider SRA -o ./input/'