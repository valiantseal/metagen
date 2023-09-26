df = read.csv('snvs_comb_res/metaseq_LogPredict_scale_covFilt_97_v2.csv')
metadat = read.csv("Final_vaxbt_dataset_AP_metadata.csv")

metadatFilt = unique(metadat[, c("AP_lab_id", "Ct_value", "WHO_variant", "days_since_last_vax", "days_post_symptom_onset", "disease_severity" )])
metadatFilt[metadatFilt == ""] <- NA
metadatFilt[metadatFilt == "."] <- NA
#metadatFilt = metadatFilt[!is.na(metadatFilt$WHO_variant),]

metadatFilt$Sample1 = gsub('_', "-", metadatFilt$AP_lab_id)

combDat = plyr::join(df, metadatFilt, by = "Sample1", type = "left", match = 'all')

# Columns to check for normality
columns_to_check <- c( "Ct_value", "days_since_last_vax", "days_post_symptom_onset", "disease_severity")

# Loop through the columns and perform Shapiro-Wilk test
for (col in columns_to_check) {
  p_value <- shapiro.test(combDat[[col]])$p.value
  
  cat(paste("Shapiro-Wilk test for column", col, "\n"))
  cat("p-value:", p_value, "\n")
  
  if (p_value > 0.05) {
    cat("Conclusion: Data is approximately normally distributed (p > 0.05)\n\n")
  } else {
    cat("Conclusion: Data is not normally distributed (p <= 0.05)\n\n")
  }
}


combDat[, columns_to_check] <- scale(combDat[, columns_to_check])

## calculate entropy

getEntropy = function(df) {
  samplesList = unique(df$Sample1)
  for (curSample in samplesList) {
    
  }
}

