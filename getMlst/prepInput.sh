#cp ../Hflu_sequences/original_files/*/*.gz ./input/
cp ../Hflu_sequences/newdir.txt ./

for i in $(cat newdir.txt)
do
mkdir -p process/"$i"
mv ./input/"$i"* ./process/"$i"/
echo "$i" > ./process/"$i"/sample.name
done
