conda activate pyseer

# similarity matrix for lmm
python ~/pyseer/scripts/phylogeny_distance.py \
--lmm /home/ubuntu/extraVol/GWAS/young2019_pyomyositis/pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/core_genome_clonMl.labelled_tree.newick \
> pyseer_input/phylogeny_kin_mat.tsv

# distance matrix for assosiation model
python ~/pyseer/scripts/phylogeny_distance.py \
/home/ubuntu/extraVol/GWAS/young2019_pyomyositis/pangenome_pirate/bactopia-runs/pangenome-20230825-133952/clonalMl/core_genome_clonMl.labelled_tree.newick \
> pyseer_input/phylogeny_dists.tsv