dir.create('./output/')

library(dplyr)
library(stringdist)

# Set the working directory to ./process/R1_001.fastq.gz
setwd("./process/R1_001.fastq.gz")

# Read the data from blastResTop_1.tsv and krakenSelVirReads.tsv
blast_df <- read.delim("blastResTop_1.tsv", stringsAsFactors = FALSE, sep = "\t")
kraken_df <- read.delim("krakenSelVirReads.tsv", stringsAsFactors = FALSE, sep = "\t")

# Merge the data frames based on the "Read" column
df <- inner_join(blast_df, kraken_df, by = "Read")

# Select and rearrange the desired columns
df <- df %>%
  select(Read, stitle, Virus)

df$stitle <- gsub("^\\s+|\\s+$|virus.*$", "virus", trimws(df$stitle), ignore.case = TRUE)
df$Virus <- gsub("^\\s+|\\s+$|virus.*$", "virus", trimws(df$Virus), ignore.case = TRUE)

# Function to calculate string similarity using Jaccard distance
calculateSimilarity <- function(stitle, Virus) {
  similarity <- stringdist::stringdistmatrix(stitle, Virus, method = "jaccard")
  return(min(similarity) < 0.8)  # Adjust the threshold as needed
}

# Apply the similarity function to each row
df$Similarity <- mapply(calculateSimilarity, df$stitle, df$Virus)

# Create a summary table with the actual virus names
summary_table <- df %>%
  group_by(Similarity, Virus) %>%
  summarize(Frequency = n()) %>%
  mutate(Virus = ifelse(Similarity, Virus, "Mismatch")) %>%
  select(Virus, Frequency)

# Consolidate all "Mismatch" entries into a single "Mismatch" row
summary_table <- summary_table %>%
  group_by(Virus) %>%
  summarize(Frequency = sum(Frequency)) %>%
  mutate(Virus = ifelse(Virus != "Mismatch", Virus, "Mismatch")) %>%
  filter(!(Virus == "Mismatch" & Frequency == 0))

# Rename and create the new "Virus" column
df <- df %>%
  rename(Blast = stitle, Kraken = Virus) %>%
  mutate(Virus = ifelse(Similarity, Kraken, "Mismatch"))

setwd("../../output")

# Save the merged data frame to a new CSV file
write.csv(summary_table, "krakBlastConfReads_summary.csv", row.names = FALSE)

# Save the modified DataFrame to a new CSV file
write.csv(df, "krakBlastConfReads.csv", row.names = FALSE)
