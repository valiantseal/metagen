conda activate bactopia

organism="Haemophilus influenzae"

bactopia --samples fastqs.txt \
	--datasets /home/ubuntu/bactopia/datasets \
	--species "$organism" \
	--coverage 150 \
	--genome_size median \
	--outdir ./bactopia_output \
	--max_cpus 4 \

