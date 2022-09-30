project="BWH_mNGS"
data_dir="2022-02-25"

mkdir input
mkdir work


conda activate dnaNexus
dx select "$project"

dx download "$data_dir"/1_raw_data/B2E* -o ./input

cd input; for i in *R1_001.fastq.gz; do echo "$i" >> ../inputList.txt; done; cd ../


cat inputList.txt |cut -f1,2,3 -d"_" > newdir.list

for i in $(cat newdir.list); do mkdir ./work/"$i"/;
cp ./input/"$i"* ./work/"$i"/;
echo "$i" > ./work/"$i"/sample.name;
done