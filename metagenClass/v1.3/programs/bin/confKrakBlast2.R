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

# Function to check if both Kraken and BLAST outputs contain either the original or the standardized virus term
check_virus <- function(blast, kraken, original_term, standardized_term) {
  if ((grepl(original_term, blast, ignore.case = TRUE) & grepl(original_term, kraken, ignore.case = TRUE)) |
      (grepl(standardized_term, blast, ignore.case = TRUE) & grepl(standardized_term, kraken, ignore.case = TRUE))) {
    return(standardized_term)
  } else {
    return(NA)
  }
}

# Apply the function for each virus term and each row of the data
for (i in 1:nrow(term_mapping)) {
  df[paste0(term_mapping$Standardized_Term[i], "_match")] <- mapply(
    check_virus, df$Blast, df$Kraken,
    MoreArgs = list(original_term = term_mapping$Original_Term[i], standardized_term = term_mapping$Standardized_Term[i])
  )
}

# Create a new column 'Virus' initialized to 'Mismatch'
df$Virus <- "Mismatch"

# Update the 'Virus' column based on matches
for (i in 1:nrow(df)) {
  for (j in 1:nrow(term_mapping)) {
    match_column <- paste0(term_mapping$Standardized_Term[j], "_match")
    if (!is.na(df[[match_column]][i]) && df[[match_column]][i] != "") {
      df$Virus[i] <- term_mapping$Standardized_Term[j]
      break
    }
  }
}

# Write the results to a new CSV
write.csv(df %>% select(Read, Blast, Kraken, Virus), "krakBlastConfReads.csv", row.names = FALSE)

# Count mismatches
num_mismatches <- nrow(df[df$Virus == "Mismatch", ])

# Gather the matches for a summary excluding mismatches
output <- df %>%
  filter(Virus != "Mismatch") %>%
  group_by(Virus) %>%
  summarize(Frequency = n()) %>%
  arrange(-Frequency)

# Append the mismatches to the summary
mismatch_df <- data.frame(Virus = "Mismatch", Frequency = num_mismatches)
output <- rbind(output, mismatch_df)

# Write the summary to a CSV
write.csv(output, "krakBlastConfReads_summary.csv", row.names = FALSE)

# Filter out mismatches
mismatches_only <- df %>%
  filter(Virus == "Mismatch")

# Write the mismatches to a CSV
write.csv(mismatches_only, "krakBlastMismatches.csv", row.names = FALSE)
