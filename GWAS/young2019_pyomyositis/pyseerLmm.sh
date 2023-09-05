mkdir -p output_gwas

conda activate pyseer

pyseer --lmm --phenotypes pyseer_input/bin_phenotype.tsv \
--kmers pyseer_input/unitig_caller.pyseer --uncompressed \
--similarity pyseer_input/phylogeny_kin_mat.tsv --output-patterns untig_patterns.txt \
--cpu 8 > ./output_gwas/lmm_unitig_pyseer.txt