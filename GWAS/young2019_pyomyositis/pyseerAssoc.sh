mkdir -p output_gwas

conda activate pyseer

pyseer --phenotypes pyseer_input/bin_phenotype.tsv \
--kmers pyseer_input/unitig_caller.pyseer --uncompressed \
--distances pyseer_input/phylogeny_dists.tsv \
--cpu 8 > ./output_gwas/assoc_unitig_pyseer.txt