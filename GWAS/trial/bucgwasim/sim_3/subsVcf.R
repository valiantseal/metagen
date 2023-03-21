library(vcfR)
vcf <- read.vcfR( '/home/ubuntu/extraVol/GWAS/bucgwassim_2/results_BacGWASim/simulations/genSim/sims.vcf', verbose = FALSE )

nrow(vcf@fix)

vcf@fix[[2,1]]

colnames(vcf@fix)

colnames(vcf@gt)
rownames(vcf@gt)

class(vcf@gt)

vcf@gt[2,3]

length(unique(vcf@gt[, 1]))

df<-extract.gt(vcf)

vcf[1:4,]

C = c(1,5,6)

vcf[C,]

z1 = c("2937", "2940" , "2957")

q1 = vcf[(vcf@fix[,2] %in% z1 )]

