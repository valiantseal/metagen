library(dplyr)
library(tidyr)
library(stringr)

# Load data
krakblast <- read.csv("output/krakBlastConfReads.csv")
library_data <- read.csv("read_label_library.csv")

# Extract main read before the '/'
krakblast$main_read <- str_extract(krakblast$Read, "^[^/]+")
library_data$main_read <- str_extract(library_data$ReadID, "^[^/]+")

# Step 1: Compare Paired-End Reads in krakBlastConfReads.csv
krakblast_grouped <- krakblast %>% 
  group_by(main_read) %>% 
  summarise(virus = toString(unique(Virus))) %>% 
  mutate(result = ifelse(str_count(virus, ",") > 0, "mismatch", virus))

# Step 2: Join the new df and read_label_library.csv
final_data <- left_join(library_data, krakblast_grouped, by = "main_read") %>% 
  select(main_read, Label, result) %>% 
  mutate(final_result = case_when(
    is.na(result) ~ "Non-Viral",
    result == Label ~ as.character(Label),
    TRUE ~ "mismatch"
  ))

# Extract only mismatched data
mismatch_data <- final_data %>% filter(final_result == "mismatch")

# Step 3: Create a summary
summary_data <- final_data %>% 
  group_by(final_result) %>% 
  summarise(count = n()) %>% 
  rename(result = final_result)

# Save the final data, mismatch data, and summary data
write.csv(final_data, "final_data.csv", row.names = FALSE)
write.csv(mismatch_data, "mismatch_data.csv", row.names = FALSE)
write.csv(summary_data, "summary_finalmatch.csv", row.names = FALSE)
