conda activate seqtk

for i in $(cat newdir.list)
do
cd process/"$i"

rm -rf splitSeq10K

mkdir splitSeq10K

seqkit split selectKraken.reads -s 900 -j 8 -O ./splitSeq10K/
cd splitSeq10K
for j in *.reads; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/reads.fas; done

cd ../../../
done
