conda activate bactopia-dev

bactopia prepare --path fastqs/ --species "Staphylococcus aureus" > fastqs.txt

bactopia --samples fastqs.txt \
--species "Staphylococcus aureus" \
--outdir ./bactopia_output \
--max_cpus 8