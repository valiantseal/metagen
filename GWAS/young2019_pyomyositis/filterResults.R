lmmDf = read.delim("output_gwas/lmm_unitig_pyseer.txt")
assocDf = read.delim("output_gwas/assoc_unitig_pyseer.txt")
lmmSnps = read.delim("output_gwas/lmm_core_pirate_vcf_pyseer.txt")
assocSnps = read.delim("output_gwas/assoc_core_pirate_vcf_pyseer.txt")


adjPFilter = function(df) {
  dfFilt = df[df$notes == "",]
  dfFilt$Bonfer_p = p.adjust(dfFilt$lrt.pvalue, method = 'bonferroni')
  dfSign = dfFilt[dfFilt$Bonfer_p < 0.05,]
  rownames(dfSign) = NULL
  return(dfSign)
}

lmmSign = adjPFilter(df = lmmDf)
assocSign = adjPFilter(df = assocDf)
lmmSnpsSign = adjPFilter(df = lmmSnps)
assocSnpsSign = adjPFilter(df = assocSnps)

overlSign = lmmSign[lmmSign$variant%in%assocSign$variant,]

lmmSignUn = lmmSign$variant
write.table(lmmSignUn, "annotate/pyseer_signUnitig.lmm", row.names = F, col.names = F, quote = F)