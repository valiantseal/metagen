conda activate bactopia

cd bactopia
for i in $(ls -d *); do cd ./"$i"
mkdir -p output
organism=$(cat bacteria.id)

bactopia --samples fastqs.txt \
	--datasets /home/ubuntu/bactopia/datasets \
	--species "$organism" \
	--coverage 150 \
	--genome_size median \
	--outdir ./output \
	--max_cpus 8 | tee bactopiaRun.log

cd ../
done
