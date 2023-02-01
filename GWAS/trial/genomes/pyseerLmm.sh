conda activate pyseer

pyseer --lmm --phenotypes ./gwas_input/pheno.tsv --continuous \
--kmers ./gwas_input/unitig_caller.pyseer --uncompressed \
--similarity gwas_input/phylogeny_kin_mat.tsv --output-patterns untig_patterns.txt \
--cpu 14 > ./gwas_output/lmm_unitig_pyseer.txt