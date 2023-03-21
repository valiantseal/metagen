conda activate pyseer

# similarity matrix for lmm
python ~/pyseer/scripts/phylogeny_distance.py --lmm results_BacGWASim/simulations/genSim/phylogeny.nwk > phylogeny_kin_mat.tsv

python ~/pyseer/scripts/phylogeny_distance.py results_BacGWASim/simulations/genSim/phylogeny.nwk > phylogeny_dists.tsv