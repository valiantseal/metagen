mkdir -p ./blastNtSummary
output_dir=$(readlink -f ./blastNtSummary)


for i in $(cat newdir.list);
do
cd process/"$i"/splitSeq10K

inst=$(ls | wc -l)

ls -d */ | parallel -j 40 'cd {} && sh ../../../../programs/bin/selVirusBlastNt0.1.sh'

cd ../../../
done