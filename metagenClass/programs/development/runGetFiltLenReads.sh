# looks like 96 cores can handle more than 60 processes
for i in $(cat newdir.list);
do
cd process/"$i"/

rm -rf virReadsFiltLen

mkdir -p virReadsFiltLen


cd splitSeq10K

ls -d */ | parallel -j 58 'cd {} && sh ~/github/DailyWork/metagenClass/programs/development/getFiltLenReads.sh'

cd ../

for j in $(cat ./lenFiltVir.reads)
do

cd splitSeq10K

for direct in *.reads_dir
do
cd "$direct"
cat virReadsFiltLen/"$j".par >> ../../virReadsFiltLen/"$j".par
cd ../
done

cd ../

done

echo "$i" "done"


cd ../../
done