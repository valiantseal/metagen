
conda activate bactopia

bactopia --samples fastqs.txt \
--datasets /home/ubuntu/extraVol/bactopia/datasets \
--species "Escherichia-coli" \
--coverage 150 \
--genome_size median \
--outdir ./bactopia_output \
--max_cpus 4 \
-profile docker -with-docker bactopia/bactopia