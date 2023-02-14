conda activate ivar

cd ./output_mpox_test/variants/bowtie2

mkdir relaxIvar

for i in $(ls *_trim.sorted.bam)
do
samtools mpileup -aa -A -d 29000000 -B -Q 0 \
"$i" | ivar variants \
-p "$i"_snps -q 15 -t 0 \
-r ~/references/idt_ref_mpox.fasta \
done

mv *_snps.tsv ./relaxIvar

cp -r relaxIvar ../../