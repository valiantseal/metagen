# Required Libraries
library(tidyverse)
library(tidytext)

# Setting the working directory and reading in the data
setwd("./process/R1_001.fastq.gz")

blast_df <- read.delim("blastResTop_1.tsv", stringsAsFactors = FALSE, sep = "\t")
kraken_df <- read.delim("krakenSelVirReads.tsv", stringsAsFactors = FALSE, sep = "\t")

setwd("../../output")

df <- inner_join(blast_df, kraken_df, by = "Read")
df <- df %>%
  select(Read, stitle, Virus) %>%
  rename(Blast = stitle, Kraken = Virus)

# Load the mapping table
term_mapping <- read.table("term_mapping.txt", stringsAsFactors = FALSE, header = TRUE)

# Function to check if the Kraken or BLAST output contains the original or the standardized virus term
find_standardized_term <- function(output, original_term, standardized_term) {
  if (grepl(original_term, output, ignore.case = TRUE) | grepl(standardized_term, output, ignore.case = TRUE)) {
    return(standardized_term)
  } else {
    return(NA)
  }
}

# Initialize KrakID and BlastID with NA
df$KrakID <- NA
df$BlastID <- NA

# Apply the function for each virus term and each row of the data to fill KrakID and BlastID
for (i in 1:nrow(term_mapping)) {
  standardized_term <- term_mapping$Standardized_Term[i]
  original_term <- term_mapping$Original_Term[i]

  # Update KrakID based on Kraken output
  krak_matches <- mapply(find_standardized_term, df$Kraken, MoreArgs = list(original_term = original_term, standardized_term = standardized_term))
  non_na_indices <- which(!is.na(krak_matches))
  df$KrakID[non_na_indices] <- krak_matches[non_na_indices]

  # Update BlastID based on Blast output
  blast_matches <- mapply(find_standardized_term, df$Blast, MoreArgs = list(original_term = original_term, standardized_term = standardized_term))
  non_na_indices <- which(!is.na(blast_matches))
  df$BlastID[non_na_indices] <- blast_matches[non_na_indices]
}

# Create a new column 'Virus' initialized to 'mismatch'
df$Virus <- "mismatch"

# Update the 'Virus' column based on matches
for (i in 1:nrow(df)) {
  if (!is.na(df$KrakID[i]) && !is.na(df$BlastID[i]) && df$KrakID[i] == df$BlastID[i]) {
    df$Virus[i] <- df$KrakID[i]
  } else {
    df$Virus[i] <- "mismatch"
  }
}

# Write the results to a new CSV
write.csv(df %>% select(Read, Blast, Kraken, KrakID, BlastID, Virus), "krakBlastConfReads.csv", row.names = FALSE)

# Count mismatches
num_mismatches <- nrow(df[df$Virus == "mismatch", ])

# Gather the matches for a summary excluding mismatches
output <- df %>%
  filter(Virus != "mismatch") %>%
  group_by(Virus) %>%
  summarize(Frequency = n()) %>%
  arrange(-Frequency)

# Append the mismatches to the summary
mismatch_df <- data.frame(Virus = "mismatch", Frequency = num_mismatches)
output <- rbind(output, mismatch_df)

# Write the summary to a CSV
write.csv(output, "krakBlastConfReads_summary.csv", row.names = FALSE)

# Filter out mismatches
mismatches_only <- df %>%

# Write the mismatches to a CSV
write.csv(mismatches_only, "krakBlastmismatches.csv", row.names = FALSE)
