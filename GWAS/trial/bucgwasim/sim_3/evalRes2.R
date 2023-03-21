adjNumb = 613431

filesList = list.files('output_gwas/lmm')

phenos = 0:9

nSignlrt<-as.numeric()
nsignBonf = as.numeric()
nsignBonf2 = as.numeric()

signlrtInRef= as.numeric()

for (i in phenos){
  inFile = paste0('./output_gwas/lmm/pheno_', i, '_pyseerVcf.txt')
  phenoFile = paste0('results_BacGWASim/simulations/phenSim/', i, '/phenSim.par')
  ref<-read.delim(phenoFile, T)
  ref$QTL<-gsub(":", "_", ref$QTL)
  
  df<- read.delim(inFile, T)
  df$BonfLrt = df$lrt.pvalue * 613431
  df$BonfLrt2 = p.adjust(df$lrt.pvalue, method = 'bonferroni')
  
  signLrt = df[(df$lrt.pvalue < 0.05),]
  signBonf = df[(df$BonfLrt< 0.05),]
  signBonf2 = df[(df$BonfLrt2< 0.05),]
  
  
  nSignlrt = c(nSignlrt, nrow(signLrt))
  nsignBonf = c(nsignBonf, nrow(signBonf))
  nsignBonf2 = c(nsignBonf2, nrow(signBonf2))
  
  signlrtInRef<- c(signlrtInRef, nrow(ref[(ref$QTL%in%signLrt$variant),]))
  
}

print(mean(nSignlrt))
print(mean(signlrtInRef))



print(mean(nsignBonf))
print(mean(nsignBonf2))

