rm -rf splitSeq

mkdir splitSeq

seqkit split virema1.fastq -s 366220 -j 8 -O ./splitSeq/

cd splitSeq

ls | wc -l 

for j in *.fastq; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/reads.fastq; done