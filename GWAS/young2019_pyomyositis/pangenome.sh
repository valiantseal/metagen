conda activate bactopia-dev

bactopia --wf pangenome --bactopia ./bactopia_output \
--outdir ./pangenome_panaroo --max_cpus 16 \
--skip_phylogeny --skip_recombination \
--run_name pirate_run \
--include pangenome_samples.txt \
--use_panaroo

