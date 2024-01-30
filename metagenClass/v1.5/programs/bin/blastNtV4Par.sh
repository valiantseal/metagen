fa_list=$(readlink -f reads.fas)
target_dir=$(pwd)
rm -f NtV4_blast.results
cd /home/ubuntu/extraVol/blast_db/nt_data/data_v4
blastn -db nt -query "$fa_list" \
-num_threads 3 -out "$target_dir"/NtV4_blast.results \
-outfmt '7 staxids qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 10
cd "$target_dir"
echo "$target_dir" "done"
