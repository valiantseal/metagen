library(doParallel)
library(foreach)

seq_info = read.csv("young2019_seqinfo.txt")
metadata = read.csv("young2019_metadata.csv")

length(metadata$Isolate.name[metadata$Isolate.name%in%seq_info$Isolate])

samplesList = unique(seq_info$Run)

fastqdl = function(curSample, provider) {
  cmd_str = paste0("fastq-dl --accession ", curSample, " --provider ", provider, " -o ./fastqs/")
  system(cmd_str)
}


dir.create("fastqs", showWarnings = F)

useCores = 6
cl <- makeCluster(useCores, type = "FORK")
registerDoParallel(cl)

combDat = foreach(i = samplesList) %dopar% {
  fastqdl(curSample = i, provider = "SRA")
}

parallel::stopCluster(cl)

