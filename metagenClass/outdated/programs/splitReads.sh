conda activate seqkit

for i in $(cat newdir.list)
do
cd process/"$i"


mkdir splitSeq10K
seqkit split krakUniq_classified.reads -s 10000 -j 8 -O ./splitSeq10K/
cd splitSeq10K
for j in *.reads; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/reads.fas; done

cd ../../../
done