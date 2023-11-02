library(ShortRead)

# List all fastq files in the 'input' directory
fastq_files <- list.files("input", pattern = "*_R1.fastq", full.names = TRUE)

# Initialize a data frame to store the read IDs and their labels
library_df <- data.frame(ReadID = character(), Label = character())

# Loop through each fastq file
for (fastq in fastq_files) {
  # Extract label from the filename
  label <- gsub("^(.*)_reads_R1\\.fastq$", "\\1", basename(fastq))

  # Read the sequences from the FASTQ
  reads <- readFastq(fastq)

  # Get the IDs of the sequences
  ids <- id(reads)

  # Combine the read IDs and labels into a data frame
  temp_df <- data.frame(ReadID = ids, Label = label)

  # Append to the master library data frame
  library_df <- rbind(library_df, temp_df)
}

# Write the library data frame to a CSV file
write.csv(library_df, "read_label_library.csv", row.names = FALSE)
