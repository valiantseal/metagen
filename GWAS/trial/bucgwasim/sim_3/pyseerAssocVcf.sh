conda activate pyseer

mkdir -p ./output_gwas

pyseer --phenotypes ./input_gwas/phnotype_0.tsv --continuous \
--vcf ./input_gwas/rename1.vcf \
--distances ./input_gwas/phylogeny_dists.tsv \
--cpu 8 > ./output_gwas/simple_assoc_pyseerVCF.txt