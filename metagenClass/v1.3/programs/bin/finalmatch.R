library(dplyr)
library(tidyr)
library(stringr)

# Load data
krakblast <- read.csv("output/krakBlastConfReads.csv")
library_data <- read.csv("read_label_library.csv")

# Extract main read before the '/'
krakblast$main_read <- str_extract(krakblast$Read, "^[^/]+")
library_data$main_read <- str_extract(library_data$ReadID, "^[^/]+")

# Compare Paired-End Reads in krakBlastConfReads.csv
krakblast_grouped <- krakblast %>%
  group_by(main_read) %>%
  summarise(
    virus = toString(unique(Virus)),
    BlastID = first(BlastID),
    KrakID = first(KrakID)
  ) %>%
  mutate(result = ifelse(str_count(virus, ",") > 0, "Mismatch", virus))

# Join the new df and read_label_library.csv
final_data <- left_join(library_data, krakblast_grouped, by = "main_read") %>%
  select(main_read, Label, result, BlastID, KrakID) %>%
  mutate(final_result = case_when(
    Label == "Non-Viral" ~ "Non-Viral",
    result == Label ~ as.character(Label),
    TRUE ~ "Mismatch"
  ))

# calculate metrics
calculate_metrics <- function(label_col, pred_col) {
  # True Positives (TP): Label matches the prediction and is not NA.
  TP <- sum(final_data[[label_col]] == final_data[[pred_col]] & !is.na(final_data[[label_col]]), na.rm = TRUE)

  # False Positives (FP): Prediction does not match Label, is not Mismatch/Non-Viral/NA, and Label is not NA.
  FP <- sum(final_data[[pred_col]] != final_data[[label_col]] &
            final_data[[pred_col]] != "Mismatch" &
            final_data[[pred_col]] != "Non-Viral" &
            !is.na(final_data[[pred_col]]) &
            !is.na(final_data[[label_col]]), na.rm = TRUE)

  # False Negatives (FN): Prediction is Mismatch/Non-Viral/NA, and Label is not NA.
  FN <- sum((final_data[[pred_col]] == "Mismatch" |
             final_data[[pred_col]] == "Non-Viral" |
             is.na(final_data[[pred_col]])) &
            !is.na(final_data[[label_col]]), na.rm = TRUE)
  
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  recall <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  f1_score <- ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)

  return(c(precision, recall, 1, f1_score))
}

# Calculate metrics for Kraken
kraken_metrics <- calculate_metrics('Label', 'KrakID')

# Calculate metrics for Blast
blast_metrics <- calculate_metrics('Label', 'BlastID')

# Calculate metrics for the Default method (which we assume is the final_result column)
default_metrics <- calculate_metrics('Label', 'final_result')

# Create a dataframe for the metrics
metrics_data <- data.frame(
  Metric = c('Precision', 'Recall', 'Specificity', 'F1'),
  Default = default_metrics,
  Kraken = kraken_metrics,
  Blast = blast_metrics,
  stringsAsFactors = FALSE
)

# Write the data frame to a CSV file
write.csv(metrics_data, "pipeline_metrics.csv", row.names = FALSE, quote = FALSE)
