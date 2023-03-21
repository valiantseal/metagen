bash -i ./programs/simmulate.sh

Rscript --vanilla ./programs/prepData.R

bash -i ./programs/unitig.sh

time bash -i ./programs/distMat.sh

time bash -i ./programs/pyseerAssoc.sh 

time bash -i ./programs/pyseerLmm.sh 

# difficult to identify a match between simulation and gwas using untigs try with SNPs

Rscript --vanilla ./programs/renameVcf.R

time bash -i ./programs/pyseerAssocVcf.sh 

time bash -i ./programs/pyseerLmmVcf.sh 