conda activate ivar

cd ./output/variants/bowtie2


ls *.sorted.bam | parallel -j 45 'samtools mpileup -aa -A -d 29000000 -B -Q 0 {} | ivar variants \
-p {}_snps -q 15 -t 0 \
-r /home/ubuntu/virilicon/references/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.fna \
-g /home/ubuntu/virilicon/references/MN908947.3/GCA_009858895.3_ASM985889v3_genomic.200409.gff'

mkdir relaxIvar
mv *_snps.tsv relaxIvar

rm -rf ../../relaxIvar
cp -r relaxIvar ../../