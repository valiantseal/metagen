
for i in $(cat newdir.list);
do
cd process/"$i"/

rm -rf virReadsFiltLen

mkdir -p virReadsFiltLen


cd splitSeq10K

time ls -d */ | parallel -j 30 'cd {} && sh ~/github/DailyWork/metagenClass/programs/development/getFiltLenReads.sh'

cd ../

for j in $(cat ./lenFiltVir.reads)
do
cat ./splitSeq10K/*/virReadsFiltLen/"$j".par >> ./virReadsFiltLen/"$j".par
done

echo "$i" "done"


cd ../../
done