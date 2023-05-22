library(biomaRt)

targDir = './Paper_figs/Fig2/'
curDate<-Sys.Date()

groupMarkers = read.csv("allDiffExprLogfc0.25_ContrVsStress_2022-10-25.csv")

# Connect to the appropriate biomart database (e.g., Ensembl)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Convert gene symbols to UniProt protein IDs
gene_symbols <- unique(groupMarkers$Genes)  # Replace with your gene symbols of interest


convert <- getBM(attributes = c("mgi_symbol", "uniprotswissprot", "uniprotsptrembl"), 
                 filters = "mgi_symbol", 
                 values = gene_symbols, 
                 mart = ensembl)


colnames(convert)[1] = "Genes"

annot = plyr::join(groupMarkers, convert, type = "left", by = "Genes", match = "all")
write.csv(annot, paste0(targDir, 'allDiffExprLogfc0.25_ContrVsStress_Uniprot_', curDate, ".csv"), row.names = F)

