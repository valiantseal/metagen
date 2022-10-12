tr ' ' '_' <virus.list > virusList.edit

mkdir -p target_results


for i in $(cat newdir.list); do rm "$i"_target.viruses; cat ./process/"$i"/splitSeq10K/*/*.par >> ./blastNtSummary/"$i"_target.viruses; done 

cd blastNtSummary

for i in $(cat ../virusList.edit)
do
for sample in $(cat ../newdir.list)
do
echo "$i" > current.vir

sed -i "s/_/ /g" current.vir
virus=$(sed -n '1p' current.vir)

grep -i "$virus" "$sample"_target.viruses > ./target_results/"$sample"__"$i"
done
done



