bash -i ./programs/simmulate.sh

Rscript --vanilla ./programs/prepData.R

time bash -i ./programs/distMat.sh

time bash -i ./programs/pyseerLmmVcf.sh

#time bash -i ./programs/pyseerEnetVcf.sh