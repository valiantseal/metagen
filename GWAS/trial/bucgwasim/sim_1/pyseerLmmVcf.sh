conda activate pyseer

pyseer --lmm --phenotypes ./input_gwas/phnotype_0.tsv --continuous \
--vcf ./input_gwas/rename1.vcf \
--similarity input_gwas/phylogeny_kin_mat.tsv \
--cpu 8 > ./output_gwas/lmm_pyseerVcf.txt