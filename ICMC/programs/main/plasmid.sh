conda activate bactopia

cd bactopia_gtdbtk
for i in $(cat ../bacteria.list); do cd ./"$i"
mkdir plasmidOut

bactopia --wf plasmidfinder \
  --bactopia ./output \
  --outdir ./plasmidOut

cd ../
done