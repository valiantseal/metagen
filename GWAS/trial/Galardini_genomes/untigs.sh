conda activate pyseer

mkdir references

ncbi-genome-download  -s refseq -A GCF_000005845.2 --assembly-levels complete --formats all -o references  --flat-output  bacteria

gzip -d ./fasta/*

(cd ./fasta && ls *) > sample.list

sed -i 's/\.fna.gz//g' sample.list

rm -f untig.input
for i in $(cat sample.list)
do
gzip -d ./bactopia_output/"$i"/assembly/"$i".fna.gz
readlink -f ./bactopia_output/"$i"/assembly/"$i".fna >> untig_input.txt
done

gzip -d references/GCF_000005845.2_ASM584v2_genomic.fna.gz
readlink -f references/GCF_000005845.2_ASM584v2_genomic.fna > untig_ref.txt

# unitig-caller is a new version of the package unitig-counter
time unitig-caller --call --refs untig_ref.txt --reads untig_input.txt --pyseer --threads 13 # 2.31 ok usage of cores

mv unitig_caller.pyseer ./gwas_input/




