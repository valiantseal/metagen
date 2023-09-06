conda activate pyseer

# make vscf file

pyseer  --phenotypes pyseer_input/bin_phenotype.tsv \
--vcf pyseer_input/pirate_core.vcf \
--distances pyseer_input/phylogeny_dists.tsv \
--cpu 8 > ./output_gwas/assoc_core_pirate_vcf_pyseer.txt