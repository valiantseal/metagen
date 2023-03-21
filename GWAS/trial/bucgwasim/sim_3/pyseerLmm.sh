conda activate pyseer

pyseer --lmm --phenotypes ./input_gwas/phnotype_0.tsv --continuous \
--kmers ./input_gwas/unitig_caller.pyseer --uncompressed \
--similarity input_gwas/phylogeny_kin_mat.tsv --output-patterns untig_patterns.txt \
--cpu 8 > ./output_gwas/lmm_unitig_pyseer.txt