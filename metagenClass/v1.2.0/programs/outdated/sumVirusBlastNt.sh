mkdir -p ./blastNtSummary
output_dir=$(readlink -f ./blastNtSummary)


for i in $(cat newdir.list);
do
#rm "$output_dir"/"$i".txt
cd process/"$i"/splitSeq10K


ls -d */ | parallel -j 12 'cd {} && sh ../../../../programs/bin/selVirusBlastNt.sh'

#for dir in $(ls -d */)
#do
#cd "$dir"
#wc -l *.sel >> "$output_dir"/"$i".txt
#cd ../
#done 

cd ../../../
done



