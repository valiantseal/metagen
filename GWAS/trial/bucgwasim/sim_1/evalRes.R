df<- read.delim('./output_gwas/lmm_unitig_pyseer.txt', T)

df$Bonfer_p<-p.adjust(df$lrt.pvalue, method = 'bonferroni')

dfFilt<-df[(df$Bonfer_p < 0.05),]
nrow(dfFilt)


dfAsoc<-read.delim('./output_gwas/simple_assoc_pyseer.txt', T)

dfAsoc$Bonfer_p<-p.adjust(dfAsoc$lrt.pvalue, method = 'bonferroni')

assocFilt<-dfAsoc[(dfAsoc$Bonfer_p < 0.05),]
nrow(assocFilt)

# vcf

df<- read.delim('./output_gwas/lmm_pyseerVcf.txt', T)

df$Bonfer_p<-p.adjust(df$lrt.pvalue, method = 'bonferroni')

dfFilt<-df[(df$Bonfer_p < 0.05),]
nrow(dfFilt)

dfAsoc<-read.delim('./output_gwas/simple_assoc_pyseerVCF.txt', T)

dfAsoc$Bonfer_p<-p.adjust(dfAsoc$lrt.pvalue, method = 'bonferroni')

assocFilt<-dfAsoc[(dfAsoc$Bonfer_p < 0.05),]
nrow(assocFilt)

# compare with original

ref<-read.delim('results_BacGWASim/simulations/phenSim/0/phenSim.par', T)

ref$QTL<-gsub(":", "_", ref$QTL)

refFilt<-ref[(ref$QTL%in%dfFilt$variant),]
nrow(refFilt)