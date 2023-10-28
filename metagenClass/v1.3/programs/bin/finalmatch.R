library(dplyr)

# Read in the data
krakBlastConfReads <- read.csv("output/krakBlastConfReads.csv", stringsAsFactors = FALSE)
read_label_library <- read.csv("read_label_library.csv", stringsAsFactors = FALSE)

# Merge the two data frames by the ReadID
merged_df <- merge(krakBlastConfReads, read_label_library, by.x = "Read", by.y = "ReadID", all.x = TRUE)

# Compare "Virus" and "Label" and create the "FinalMatch" column
merged_df$FinalMatch <- ifelse(merged_df$Virus == merged_df$Label, merged_df$Virus, "Mismatch")

# Create a summary dataframe with the counts of each "FinalMatch"
summary_df <- merged_df %>%
  group_by(FinalMatch) %>%
  summarize(Count = n())

# Save the summary to a CSV file
write.csv(summary_df, "summary_finalmatch.csv", row.names = FALSE)
