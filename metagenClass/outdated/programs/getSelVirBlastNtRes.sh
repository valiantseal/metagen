mkdir -p krUnVipr/blastNtSummary/results
output_dir=$(readlink -f krUnVipr/blastNtSummary/results)

for i in $(cat newdir.list);
do
cat krUnVipr/"$i"/all_classified_reads/splitSeq10K/*/alphaherpesvirus.sel >> "$output_dir"/"$i"_alphaherpesvirus.txt
cat krUnVipr/"$i"/all_classified_reads/splitSeq10K/*/mastadenovirus.sel >> "$output_dir"/"$i"_mastadenovirus.txt
cat krUnVipr/"$i"/all_classified_reads/splitSeq10K/*/polyomavirus.sel >> "$output_dir"/"$i"_polyomavirus.txt
done

