combStats = read.csv("ampseq_metaseq_overlap_combStats.csv")

dir.create("Overlap_Pos_mapping/Comb_fasta/", recursive = T, showWarnings = F)

makeAlignment = function(combStats, i) {
  ampSample = combStats$OrigName[i]
  metaSample = combStats$OrigName_meta[i]
  sampName = combStats$Sample[i]
  catFile = paste0("Overlap_Pos_mapping/Comb_fasta/",   sampName, ".fasta")
  ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
  metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
  cmd_str = paste("cat references/MN908947.3.fna", ampPath, metaPath, ">", catFile, sep = " ")
  system(cmd_str)
  cmd_mafft = paste0("mafft --maxiterate 1000 --thread 4 --globalpair ",  catFile, " > Overlap_Pos_mapping/", sampName, "_align.fasta")
  system(cmd_mafft)
}

makeAlignment(combStats=combStats, i=1)

dir.create("Overlap_Pos_mapping/Indecies/", recursive = T, showWarnings = F)


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


compileIndex = function(combStats, i) {
  ampSample = combStats$OrigName[i]
  metaPath = combStats$OrigName_meta[i]
  sampName = combStats$Sample[i]
  ampPath = paste0("IntraSnv_ampseq_overlap/",  ampSample, "/reference.fa")
  metaPath = paste0("IntraSnv_metaseq_overlap/",  metaSample, "/reference.fa")
  alignmentPath = paste0("Overlap_Pos_mapping/", sampName, "_align.fasta")
  
  # read fasta
  curAlignment = phylotools::read.fasta(alignmentPath)
  ampFasta = phylotools::read.fasta(ampPath)
  metaFasta = phylotools::read.fasta(metaPath)
  
  # split strings to vectors
  ampAlignVect = toupper(strsplit(curAlignment[2,2], "")[[1]])
  metaAlignVect = toupper(strsplit(curAlignment[3,2], "")[[1]])
  
}