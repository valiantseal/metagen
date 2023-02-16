rm -rf splitSeq

mkdir splitSeq

seqkit split virema1.fastq -s 1000000 -j 4 -O ./splitSeq/

cd splitSeq

for j in *.fastq; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/reads.fastq; done