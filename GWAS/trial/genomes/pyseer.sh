conda activate pyseer

mkdir -p ./gwas_output

pyseer --phenotypes ./gwas_input/pheno.tsv --continuous \
--kmers ./gwas_input/unitig_caller.pyseer --uncompressed \
--distances ./gwas_input/phylogeny_dists.tsv \
--cpu 14 > ./gwas_output/simple_assoc_pyseer.txt