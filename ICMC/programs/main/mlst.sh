conda activate bactopia

cd bactopia_gtdbtk
for i in $(cat ../bacteria.list); do cd ./"$i"
mkdir mlstOut

bactopia --wf mlst \
--bactopia ./output \
--outdir ./mlstOut \

cd ../
done