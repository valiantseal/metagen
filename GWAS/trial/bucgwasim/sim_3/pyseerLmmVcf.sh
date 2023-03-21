conda activate pyseer

mkdir -p output_gwas/lmm

for i in $(cd phenotypes &&  ls)
do

outF=$(echo "$i" | cut -f1 -d ".")"_pyseerVcf.txt"

pyseer --lmm --phenotypes ./phenotypes/"$i" --continuous \
--vcf /home/ubuntu/extraVol/GWAS/bucgwassim_2/results_BacGWASim/simulations/genSim/sims.vcf \
--similarity phylogeny_kin_mat.tsv \
--cpu 15 > ./output_gwas/lmm/"$outF"

done

