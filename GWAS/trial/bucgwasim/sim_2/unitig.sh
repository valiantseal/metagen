conda activate pyseer

mkdir input_gwas


unitig-caller --call --refs untig_ref.txt --reads untig_input.txt --pyseer --threads 13 # 2.31 ok usage of cores

mv unitig_caller.pyseer ./input_gwas/