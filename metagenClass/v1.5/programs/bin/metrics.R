library(dplyr)
library(ggplot2)

# Read the CSV files
results <- read.csv("fullresults_with_labels.csv", sep = "\t", header = TRUE)

# Function to calculate precision
calculate_precision <- function(df, tool_col, label_col) {
  TP <- sum(df[[tool_col]] == df[[label_col]])
  FP <- sum(df[[tool_col]] != df[[label_col]] & !is.na(df[[tool_col]]))
  
  Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  
  return(Precision)
}

# Calculate precision for each label and tool
results <- results %>%
  mutate(
    Precision_Pipeline = ifelse(LCA == Label, 1, 0),
    Precision_BLAST1_1 = ifelse(BID1_1 == Label, 1, 0),
    Precision_BLAST1_2 = ifelse(BID1_2 == Label, 1, 0),
    Precision_Kraken1_1 = ifelse(KID1_1 == Label, 1, 0),
    Precision_Kraken1_2 = ifelse(KID1_2 == Label, 1, 0)
  )

# Aggregate by Label
precision_by_label <- results %>%
  group_by(Label) %>%
  summarise(
    Precision_Pipeline = mean(Precision_Pipeline, na.rm = TRUE),
    Precision_BLAST = mean(c(Precision_BLAST1_1, Precision_BLAST1_2), na.rm = TRUE),
    Precision_Kraken = mean(c(Precision_Kraken1_1, Precision_Kraken1_2), na.rm = TRUE)
  )

# Plot Precision for each Label (Figure A)
precision_plot <- ggplot(precision_by_label, aes(x = Label)) +
  geom_bar(aes(y = Precision_Pipeline, fill = "Pipeline"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = Precision_BLAST, fill = "BLAST"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = Precision_Kraken, fill = "Kraken"), stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Pipeline" = "blue", "BLAST" = "red", "Kraken" = "green")) +
  theme_minimal() +
  labs(y = "Precision", title = "Precision by Label for Each Tool")

# Calculate overall precision with standard deviation for error bars (Figure B)
overall_precision <- results %>%
  summarise(
    Precision_Pipeline = mean(Precision_Pipeline, na.rm = TRUE),
    Precision_BLAST = mean(c(Precision_BLAST1_1, Precision_BLAST1_2), na.rm = TRUE),
    Precision_Kraken = mean(c(Precision_Kraken1_1, Precision_Kraken1_2), na.rm = TRUE),
    SD_Pipeline = sd(Precision_Pipeline, na.rm = TRUE),
    SD_BLAST = sd(c(Precision_BLAST1_1, Precision_BLAST1_2), na.rm = TRUE),
    SD_Kraken = sd(c(Precision_Kraken1_1, Precision_Kraken1_2), na.rm = TRUE)
  ) %>%
  gather(key = "Tool", value = "Precision", -SD_Pipeline, -SD_BLAST, -SD_Kraken) %>%
  mutate(Error = ifelse(Tool == "Precision_Pipeline", SD_Pipeline,
                        ifelse(Tool == "Precision_BLAST", SD_BLAST, SD_Kraken)))

# Plot Overall Precision with error bars (Figure B)
overall_precision_plot <- ggplot(overall_precision, aes(x = Tool, y = Precision, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Precision - Error, ymax = Precision + Error), width = 0.25) +
  scale_fill_manual(values = c("Pipeline" = "blue", "BLAST" = "red", "Kraken" = "green")) +
  theme_minimal() +
  labs(y = "Precision", title = "Overall Precision for Each Tool")

# Save the plots
ggsave("precision_by_label_plot.png", precision_plot, width = 10, height = 6)
ggsave("overall_precision_plot.png", overall_precision_plot, width = 10, height = 6)
