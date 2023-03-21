conda activate pyseer


pyseer --wg enet --phenotypes pheno_0.tsv --continuous \
--vcf /home/ubuntu/extraVol/GWAS/bucgwassim_2/results_BacGWASim/simulations/genSim/sims.vcf \
--distances phylogeny_dists.tsv \
--cpu 14 > ./output_gwas/enet_pyseerVcf.txt