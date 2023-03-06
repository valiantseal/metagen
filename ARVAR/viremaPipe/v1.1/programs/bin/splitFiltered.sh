rm -rf splitSeq

mkdir splitSeq

lines=$(cat virema1.fastq | wc -l)

parts=$(echo "$lines"/4/95+5|bc)
echo "$parts"

seqkit split virema1.fastq -s "$parts" -j 4 -O ./splitSeq/

cd splitSeq

for j in *.fastq; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/reads.fastq; done