for sample in $(cat newdir.list);
do
for id in $(cat ./kraqSummary/kraken.ids)
do
#awk -F '\t' '{ if ($3 == "$id") { print } }' ./process/"$sample"/krakUniq_sample.kraken >> ./kraqSummary/"$sample".reads

grep "$id"  ./process/"$sample"/krakUniq_sample.kraken >> ./kraqSummary/"$sample".reads

#awk -v pat="$id" -F '\t' '$3~pat' ./process/"$sample"/krakUniq_sample.kraken >> ./kraqSummary/"$sample".reads

done
done
