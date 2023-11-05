# Required Libraries
library(tidyverse)
library(tidytext)

# Setting the working directory and reading in the data
dir.create('./output/')
setwd("./process/R1_001.fastq.gz")

blast_df <- read.delim("blastResTop_1.tsv", stringsAsFactors = FALSE, sep = "\t")
kraken_df <- read.delim("krakenSelVirReads.tsv", stringsAsFactors = FALSE, sep = "\t")
setwd("../../output")

df <- inner_join(blast_df, kraken_df, by = "Read")
df <- df %>%
  select(Read, stitle, Virus) %>%
  rename(Blast = stitle, Kraken = Virus)

# Extract terms from Blast and Kraken columns separately
blast_terms <- df %>%
  mutate(Blast_Terms = str_extract_all(Blast, boundary("word"))) %>%
  unnest(Blast_Terms) %>%
  select(Blast_Terms)

kraken_terms <- df %>%
  mutate(Kraken_Terms = str_extract_all(Kraken, boundary("word"))) %>%
  unnest(Kraken_Terms) %>%
  select(Kraken_Terms)

# Find common terms and their frequencies
common_terms <- inner_join(blast_terms, kraken_terms, by = c("Blast_Terms" = "Kraken_Terms")) %>%
  group_by(Term = Blast_Terms) %>%
  tally(name = "Frequency") %>%
  arrange(-Frequency)

# Filter out terms that are shorter than 3 characters and non-viral terms
non_viral_terms <- c("genome", "cds", "strain", "viral", "complete", "segment", "protein", "gene", "isolate", "polymerase", "syndrome", "virus", "human", "sapien", "sequence")
common_terms <- common_terms %>%
  filter(nchar(Term) > 3 & !Term %in% non_viral_terms) %>%
  head(30)

# Generate lookup table
lookup_table <- common_terms %>%
  rename(Original_Term = Term) %>%
  mutate(Standardized_Term = Original_Term)

# Save the lookup table for manual editing
write.table(lookup_table, "term_mapping.txt", row.names = FALSE, sep = "\t", quote = FALSE)

