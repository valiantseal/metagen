for i in $(cat newdir.list)
do
mkdir -p ./output/"$i"/merged_reads
mkdir -p ./output/"$i"/R1_not_merged
mkdir -p ./output/"$i"/R2_not_merged
for j in $(cat ./work/"$i"/merged.list)
do 
cp ./work/"$i"/"$j"/classified.reads ./output/"$i"/"$j"/
cp ./work/"$i"/"$j"/blast.results ./output/"$i"/"$j"/
done 
done


