mkdir -p ./blastNtSummary/results
output_dir=$(readlink -f ./blastNtSummary/results)

for i in $(cat newdir.list);
do

rm "$output_dir"/"$i"_alphaherpesvirus.txt
rm "$output_dir"/"$i"_mastadenovirus.txt
rm "$output_dir"/"$i"_polyomavirus.txt

cat process/"$i"/splitSeq10K/*/alphaherpesvirus.sel >> "$output_dir"/"$i"_alphaherpesvirus.txt
cat process/"$i"/splitSeq10K/*/mastadenovirus.sel >> "$output_dir"/"$i"_mastadenovirus.txt
cat process/"$i"/splitSeq10K/*/polyomavirus.sel >> "$output_dir"/"$i"_polyomavirus.txt
done

