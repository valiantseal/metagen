seq_info = read.csv("young2019_seqinfo.txt")
metadata = read.csv("young2019_metadata.csv")
bactRun = read.delim("fastqs.txt")

colnames(seq_info)

seqInfoSel = seq_info[, c("Run", "Isolate", "Organism")]
colnames(metadata)
colnames(metadata)[1] = "Isolate"

combMetadat = plyr::join(metadata, seqInfoSel, by = "Isolate", type = "left", match = "all")

curMetadata = combMetadat[combMetadat$Run%in%bactRun$sample,]

unique(curMetadata$Organism)

samples = curMetadata$Run

write.table(samples, "pangenome_samples.txt", row.names = F, col.names = F, quote = F)

write.csv(curMetadata, "metadata.csv", row.names = F)

binPhenotype = function(df) {
  df$Phenotype = NA
  for ( i in 1:nrow(df) ) {
    if (df$Case.control.status[i] == "Control") {
      df$Phenotype[i] = 0
    } else if (df$Case.control.status[i] == "Case") {
      df$Phenotype[i] = 1
    }
  }
  return(df)
}

curMetadata = binPhenotype(df =curMetadata)

curPheno = curMetadata[, c("Run", "Phenotype")]

write.table(curPheno, "pyseer_input/bin_phenotype.tsv", row.names = F, col.names = T, quote = F, sep = "\t")