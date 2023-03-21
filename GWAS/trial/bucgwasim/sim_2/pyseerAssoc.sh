conda activate pyseer

mkdir -p ./output_gwas

pyseer --phenotypes ./input_gwas/phnotype_0.tsv --continuous \
--kmers ./input_gwas/unitig_caller.pyseer --uncompressed \
--distances ./input_gwas/phylogeny_dists.tsv \
--cpu 8 > ./output_gwas/simple_assoc_pyseer.txt