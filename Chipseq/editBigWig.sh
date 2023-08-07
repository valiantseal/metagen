

# /home/flyhunter/Bio_bins/bigWigToBedGraph c1_scaleRc.bw c1_scaleRc.bg
# awk -F'\t' -v OFS='\t' '$1 == "chr4" { print }' c1_scaleRc.bg > c1_scaleRc_chr4.bg
# awk -F'\t' -v OFS='\t' '{ if ($4 < 0) $4 = 0; print }' c1_scaleRc_chr4.bg > c1_scaleRc_chr4_edit.bg
# /home/flyhunter/Bio_bins/bedSort c1_scaleRc_chr4_edit.bg c1_scaleRc_chr4_edit_sort.bg
# /home/flyhunter/Bio_bins/bedGraphToBigWig c1_scaleRc_chr4_edit_sort.bg ../mm10.chromsize c1_scaleRc_chr4_edit_sort.bw


for i in *scaleRc.bw;
do 
baseName="${i%.bw}"
/home/flyhunter/Bio_bins/bigWigToBedGraph "$i" inter_files/"$baseName".bg
awk -F'\t' -v OFS='\t' '$1 == "chr4" { print }' inter_files/"$baseName".bg > inter_files/"$baseName"_chr4.bg
awk -F'\t' -v OFS='\t' '{ if ($4 < 0) $4 = 0; print }'  inter_files/"$baseName"_chr4.bg >  inter_files/"$baseName"_chr4_edit.bg
/home/flyhunter/Bio_bins/bedSort inter_files/"$baseName"_chr4_edit.bg inter_files/"$baseName"_chr4_edit_sort.bg
/home/flyhunter/Bio_bins/bedGraphToBigWig inter_files/"$baseName"_chr4_edit_sort.bg ../mm10.chromsize "$baseName"_chr4_edit_sort.bw
done