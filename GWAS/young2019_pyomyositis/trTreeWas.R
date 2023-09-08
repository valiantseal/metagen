library(treeWAS)

dna <- read.dna(file = "/home/ubuntu/extraVol/GWAS/young2019_pyomyositis/pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/core-genome.aln", format = "fasta")

mat <- DNAbin2genind(dna)@tab

tree <- read.tree(file = "/home/ubuntu/extraVol/GWAS/young2019_pyomyositis/pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/core_genome_clonMl.labelled_tree.newick")


phen_df = read.csv("metadata.csv")

phen <- phen_df$Case.control.status
names(phen) <- phen_df$Run

## Cross-check labels with each other:
all(tree$tip.label %in% rownames(mat))
all(rownames(mat) %in% tree$tip.label)
all(tree$tip.label %in% names(phen))
all(names(phen) %in% tree$tip.label)
all(names(phen) %in% rownames(mat))
all(rownames(mat) %in% names(phen))

suffixes <- keepLastN(colnames(mat), n = 2)
suffixes <- unique(suffixes)

if(all(suffixes %in% c(".a", ".c", ".g", ".t"))){
  ## SNPs:
  #snps <- get.binary.snps(mat)
  snps = mat
}

#hist(phen)

dir.create("treeWAS", showWarnings = F)

out <- treeWAS(snps = snps,
               phen = phen,
               tree = tree,
               seed = 1, filename.plot = "treeWAS/tr_output_allCol.pdf")


saveRDS(out, "treeWAS/tr_output_allCol.rds")

print(out, sort.by.p=FALSE)
