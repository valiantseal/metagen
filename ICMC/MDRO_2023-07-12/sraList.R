df = read.csv("metadata/SRA_IDS_GAMuGSI_updated07052.csv")

length(unique(df$SRA.Number))

sra = unique(df$SRA.Number)

write.table(sra, "sra_samples.list", row.names = F, col.names = F, quote = F)