conda activate pyseer

# regular genetic distance
python ~/pyseer/scripts/phylogeny_distance.py ./tree_rename.nwk > input_gwas/phylogeny_dists.tsv

# similarity matrix for lmm
python ~/pyseer/scripts/phylogeny_distance.py --lmm ./tree_rename.nwk > input_gwas/phylogeny_kin_mat.tsv

