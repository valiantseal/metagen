library(treeWAS)


runTreeWas = function(fastaPath, treePath, runName) {
  dna <- read.dna(file = fastaPath, format = "fasta")
  tree <- read.tree(file = treePath)
  phen_df = read.csv("metadata.csv")
  
  mat <- DNAbin2genind(dna)@tab
  phen <- phen_df$Case.control.status
  names(phen) <- phen_df$Run
  
  suffixes <- keepLastN(colnames(mat), n = 2)
  suffixes <- unique(suffixes)
  
  if(all(suffixes %in% c(".a", ".c", ".g", ".t"))){
    ## SNPs:
    snps <- get.binary.snps(mat)
  }
  
  out <- treeWAS(snps = snps,
                 phen = phen,
                 tree = tree,
                 seed = 1, filename.plot = paste0("treeWAS/plots_", runName, ".pdf"))
  
  saveRDS(out, paste0("treeWAS/", runName, ".rds"))
  print(out, sort.by.p=FALSE)
  
}

runTreeWas(fastaPath = "pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/pangenome_alignment.fasta", 
           treePath = "pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/pangenome_clonMl.labelled_tree.newick", 
           runName = "pirate_pangenome")