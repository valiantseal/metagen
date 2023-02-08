conda activate dnaNexus

project="Dengue_virus"
dir='2022-12-08 Nica_DENV B1r'

mkdir input
mkdir process

dx select "$project"


for i in $(cat downloadSamples.list)
do
dx download "$dir"/1_raw_data/"$i"* -o ./input
done
