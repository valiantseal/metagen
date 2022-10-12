for j in $(cat merged.list)
do 
cd "$j"


sample=$(cat merge.type)

fa_list=$(readlink -f krakUniq_classified.reads)
target_dir=$(pwd)

cd /home/ubuntu/vipr

blastn -db NONFLU_All.nt -query "$fa_list" \
-num_threads 8 -out "$target_dir"/krakenUniq_Vipr_blast.results \
-outfmt '7 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore'

cd "$target_dir"

cd ../
done