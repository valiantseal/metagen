conda activate pyseer

mkdir output_gwas

pyseer --lmm --phenotypes pheno_0.tsv --continuous \
--vcf /home/ubuntu/extraVol/GWAS/bucgwassim_2/results_BacGWASim/simulations/genSim/sims.vcf \
--similarity phylogeny_kin_mat.tsv \
--cpu 14 > ./output_gwas/lmm_pyseerVcf.txt

pyseer --lmm --phenotypes pheno_0.tsv --continuous \
--vcf /home/ubuntu/extraVol/GWAS/bucgwassim_2/results_BacGWASim/simulations/genSim/sims_no_selection.vcf \
--similarity phylogeny_kin_mat.tsv \
--cpu 14 > ./output_gwas/lmm_pyseerVcfNoSel.txt