library(foreach)
library(doParallel)

combStats = read.csv("ampseq_metaseq_overlap_combStats.csv")

dir.create("Overlap_Pos_mapping/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats) {
  system("cat references/MN908947.3.fna > Overlap_Pos_mapping/Comb_fasta/allSequences.fasta")
  seqType = character()
  seqNames = character()
  for ( i in 1:nrow(combStats) ) {
    ampSample = combStats$OrigName[i]
    metaSample = combStats$OrigName_meta[i]
    sampName = combStats$Sample[i]
    seqNames= c(seqNames, ampSample, metaSample)
    seqType = c(seqType, "ampseq", "metaseq")
    
    catFile = "Overlap_Pos_mapping/Comb_fasta/allSequences.fasta"
    ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
    metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
    
    cmd_str = paste("cat ", ampPath, metaPath, ">>", catFile, sep = " ")
    system(cmd_str)
  }
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 29 --globalpair ",  catFile, " > Overlap_Pos_mapping/All_samples_align.fasta")
  system(cmd_mafft)
  
  combDat = data.frame(seqNames, seqType)
  write.csv(combDat, "All_overlap_samples_alignment_order.csv", row.names = F)
}


makeAlignment(combStats)
