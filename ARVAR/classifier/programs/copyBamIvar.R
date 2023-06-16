setwd("/home/ubuntu/extraVol/ARVAR/classifier/Vivacilty_v1.0.1")

df = read.csv("../Ludy_metaAmpIvar_overlapSnv.csv")

samplesList = unique(df$FullSamp)

#filesList = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-16/output/variants/bowtie2')

for ( sample in samplesList) {
  targDir = paste0("process_par/", sample, "/")
  dir.create(targDir, recursive = T, showWarnings = F)
  sampBam = list.files('/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-16/output/variants/bowtie2', pattern = sample, full.names = T)
  for (curBam in sampBam) {
    if (grepl(".bai", curBam)) {
      file.copy(curBam, paste0(targDir, "output.bam.bai"))
    } else {
      file.copy(curBam, paste0(targDir, "output.bam"))
    }
  
  }
}