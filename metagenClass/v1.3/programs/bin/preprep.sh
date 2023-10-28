mkdir -p input

# Loop through each fasta file and generate reads using iss
for fasta in *.fasta; do
    # Get the base name of the fasta file (without extension)
    base_name=$(basename "$fasta" .fasta)

    # Generate reads with iss
    iss generate -g "$fasta" -o "input/${base_name}_reads" -n 1000 --cpus 16 --model miseq
done

# Create the library CSV using the R script
Rscript --vanilla ./programs/bin/create_library.R

# Concatenate all _R1 reads into one combined file
cat input/*_R1.fastq > input/R1_001.fastq

# Concatenate all _R2 reads into one combined file
cat input/*_R2.fastq > input/R2_001.fastq

# Gzip both files
gzip input/R1_001.fastq
gzip input/R2_001.fastq
