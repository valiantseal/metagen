mkdir -p ./blastNtSummary
output_dir=$(readlink -f ./blastNtSummary)

for i in $(cat newdir.list);
do
cd process/"$i"/splitSeq10K

inst=$(ls | wc -l)

ls -d */ | parallel -j 95 'cd {} && Rscript --vanilla ~/github/DailyWork/metagenClass/programs/development/filterBlast.R'

cd ../../../
done