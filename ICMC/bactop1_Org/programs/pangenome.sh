conda activate bactopia


mkdir ./panOut

bactopia --wf pangenome \
--bactopia ./bactopia_output \
--outdir ./panOut \
--max_cpus 32 \
--skip_phylogeny \
--skip_recombination