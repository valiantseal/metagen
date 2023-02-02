cp /home/ubuntu/Bacgwasim/results_BacGWASim/simulations/genSim/sims.vcf ./input/snps.vcf

Rscript --vanilla ./programs/prepData.R

pyseer --phenotypes ./input/pheno.tsv --continuous \
--vcf ./input/snps.vcf \
--distances ./input/distance.tsv \
--min-af 0.01 --max-af 0.99 --cpu 8 > pyseer.assoc