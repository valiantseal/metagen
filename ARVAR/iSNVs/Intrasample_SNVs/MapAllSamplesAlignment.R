library(foreach)
library(doParallel)
library(phylotools)

combStats = read.csv("ampseq_metaseq_overlap_combStats.csv")

dir.create("Overlap_Pos_mapping/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats) {
  #system("cat references/MN908947.3.fna > Overlap_Pos_mapping/Comb_fasta/allSequences.fasta")
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
    #system(cmd_str)
  }
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 29 --auto ",  catFile, " > Overlap_Pos_mapping/All_samples_align.fasta")
  system(cmd_mafft)
  
  combDat = data.frame(seqNames, seqType)
  write.csv(combDat, "All_overlap_samples_alignment_order.csv", row.names = F)
}

# run alignment function
makeAlignment(combStats)

# make index
makeIndex = function(alignTargVec) {
  Alignment_Pos = numeric()
  Original_Pos = numeric()
  alignBases = character()
  curOrigIndex = 0
  for ( i in 1:length(alignTargVec) ) {
    curBase = alignTargVec[i]
    if ( curBase != "-" ) {
      Alignment_Pos = c(Alignment_Pos, i)
      curOrigIndex = curOrigIndex + 1
      Original_Pos = c(Original_Pos, curOrigIndex)
      #alignBases = c(alignBases,  curBase)
    }
  }
  combIndex = data.frame(Original_Pos, Alignment_Pos)
  return(combIndex)
}

# 

alignOrder = read.csv("All_overlap_samples_alignment_order.csv")
alignFasta = read.fasta("Overlap_Pos_mapping/All_samples_align.fasta")
alignFasta = alignFasta[2:nrow(alignFasta),]
alignFasta$seq.name = gsub(" MN908947.3", "", alignFasta$seq.name)
identical(alignFasta$seq.name, alignOrder$seqNames)
alignFasta$SeqType = alignOrder$seqType

compileIndex = function(combStats, i, makeCheck, alignFasta) {
  
}
