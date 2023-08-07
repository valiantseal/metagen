genes_list = read.table("/home/flyhunter/Kai/Enrichment/genes.list", T)
genes_list = toupper(genes_list$Genes)

motifs_df = read.delim("/home/flyhunter/Kai/Chipseq/pnas/process_custom/motifs/chd_all_shun/knownResults.txt")

motifs_sign = motifs_df[(motifs_df$P.value < 0.05),]

allMotifs = toupper(motifs_sign$Motif.Name)

findOverlapGenes = function(genes, motifs) {
  combGenes = character()
  for (gene in genes) {
    if (any(grepl(gene, motifs, ignore.case = T))) {
      combGenes = c(combGenes, gene)
    }
  }
  return(combGenes)
}

overlapGenes = findOverlapGenes(genes = genes_list, motifs = allMotifs)