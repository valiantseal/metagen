mkdir -p krUnVipr/blastNtSummary
output_dir=$(readlink -f krUnVipr/blastNtSummary)


for i in $(cat newdir.list);
do
cd krUnVipr/"$i"/all_classified_reads/splitSeq10K


ls -d */ | parallel -j 6 'cd {} && sh ../../../../../programs/selVirusBlastNt.sh'

for dir in $(ls -d */)
do
cd "$dir"
wc -l *.sel >> "$output_dir"/"$i".txt
cd ../
done 

cd ../../../../
done



