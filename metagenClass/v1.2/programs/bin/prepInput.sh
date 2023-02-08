
cp ./programs/virus.list ./

cd input; for i in *R1_001.fastq.gz; do echo "$i" >> ../inputList.txt; done; cd ../


cat inputList.txt |cut -f1,2,3 -d"_" > newdir.list


for i in $(cat newdir.list); do mkdir -p ./process/"$i"/;
cp ./input/"$i"* ./process/"$i"/;
echo "$i" > ./process/"$i"/sample.name;
done

