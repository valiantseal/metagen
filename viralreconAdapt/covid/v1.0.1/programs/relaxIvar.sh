conda activate ivar

cd ./output/variants/bowtie2

rm -rf relaxIvar
mkdir -p relaxIvar

for i in $(ls *_trim.sorted.bam)
do
samtools mpileup -aa -A -d 29000000 -B -Q 0 \
"$i" | ivar variants \
-p "$i"_snps -q 15 -t 0 \
-r /home/ubuntu/virilicon/references/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.fna \
-g /home/ubuntu/virilicon/references/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff
done

mv *_snps.tsv ./relaxIvar

rm -rf ../../relaxIvar
cp -r relaxIvar ../../
