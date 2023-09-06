conda activate pyseer

# make vscf file
snp-sites -v -o pyseer_input/pirate_core.vcf /home/ubuntu/extraVol/GWAS/young2019_pyomyositis/pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/core-genome.aln

pyseer --lmm --phenotypes pyseer_input/bin_phenotype.tsv \
--vcf pyseer_input/pirate_core.vcf \
--similarity pyseer_input/phylogeny_kin_mat.tsv --output-patterns vcf_core_pirate_patterns.txt \
--cpu 8 > ./output_gwas/lmm_core_pirate_vcf_pyseer.txt