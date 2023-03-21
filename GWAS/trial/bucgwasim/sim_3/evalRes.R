
# vcf

df<- read.delim('./output_gwas/lmm_pyseerVcf.txt', T)

df$Bonfer_p<-p.adjust(df$lrt.pvalue, method = 'bonferroni')

dfFilt<-df[(df$lrt.pvalue < 0.05),]


dfAsoc<-read.delim('./output_gwas/lmm_pyseerVcfNoSel.txt', T)

dfAsoc$Bonfer_p<-p.adjust(dfAsoc$lrt.pvalue, method = 'bonferroni')

assocFilt<-dfAsoc[(dfAsoc$lrt.pvalue < 0.05),]

# compare with original

ref<-read.delim('results_BacGWASim/simulations/phenSim/0/phenSim.par', T)

ref$QTL<-gsub(":", "_", ref$QTL)

refFilt<-ref[(ref$QTL%in%dfFilt$variant),]

# manual p-value adjustment
df<- read.delim('./output_gwas/lmm_pyseerVcf.txt', T)
df$Adj_p<-df$lrt.pvalue*30000
dfFilt<-df[(df$Adj_p < 0.05),]



# 
refFiltSel<-ref[(ref$QTL%in%assocFilt$variant),]

# elastci net
enet<- read.delim('./output_gwas/enet_pyseerVcf.txt', T)

enet$Bonfer_p<-p.adjust(enet$filter.pvalue, method = 'bonferroni')

enetFilt<-enet[(enet$filter.pvalue < 0.05),]

enetBonf<- enet[(enet$Bonfer_p < 0.05),]

nrow(ref[(ref$QTL%in%enetFilt$variant),])

nrow(ref[(ref$QTL%in%enetBonf$variant),])