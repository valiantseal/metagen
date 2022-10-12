mkdir -p ./blastNtSummary/results
output_dir=$(readlink -f ./blastNtSummary/results)

rm "$output_dir"/*.txt

for i in $(cat newdir.list);
do
for virus in $(cat virus.list)
do
cat process/"$i"/splitSeq10K/*/"$virus".sel >> "$output_dir"/"$i"_"$virus".txt
done
done

