# iqtree will replace - with _ so need to plan in advance if it is better to rename samples

# read.table() would replace - with . in column names (matricies), read_delim() would not

# names in the distance and kinship matricies influence results, better to rename files if use iqtree in advance

bash -i ./programs/bactopia.sh

bash -i ./programs/pangenome.sh

Rscript --vanilla ./programs/iqtree.R

Rscript --vanilla ./programs/relabelTree.R

# need untigs, distance from phylogeny, and phenotypes

Rscript --vanilla ./programs/prepData.R

mkdir gwas_input; cp pheno.tsv ./gwas_input

bash -i ./programs/untigs.sh

python ~/pyseer/scripts/phylogeny_distance.py ./iqtree/relab_core.contree > gwas_input/phylogeny_dists.tsv

time bash -i ./programs/pyseer.sh # 3.13 14cpu

python ~/pyseer/scripts/phylogeny_distance.py --lmm ./iqtree/relab_core.contree > gwas_input/phylogeny_kin_mat.tsv

time bash -i ./programs/pyseerLmm.sh

# get snps vcf file
snp-sites -v -o ./gwas_input/core_aln_smps.vcf ./panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz

snp-sites -p -o ./gwas_input/core_aln_smps.phy ./panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz