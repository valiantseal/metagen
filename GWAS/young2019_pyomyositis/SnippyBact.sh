conda activate bactopia-dev

bactopia --wf snippy --bactopia ./bactopia_output \
--outdir ./snippy_output --max_cpus 31 \
--skip_phylogeny \
--run_name "bact_snippy" \
--reference s_aureus.gb