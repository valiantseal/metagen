for j in $(cat merged.list)
do 
cd "$j"

sample=$(cat merge.type)

fa_list=$(readlink -f readsForBlast.fa)
target_dir=$(pwd)

cd /home/ubuntu/vipr

blastn -db NONFLU_All.nt -query "$fa_list" \
-num_threads 7 -out "$target_dir"/blast.results \
-outfmt '7 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-max_target_seqs 10

cd "$target_dir"

cd ../
done