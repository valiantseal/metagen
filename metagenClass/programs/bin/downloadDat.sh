conda activate dnaNexus

project="BWH_mNGS"

mkdir input
mkdir process

dx select "$project"


for i in $(cat downloadSamples.list)
do
dx download "$i"* -o ./input
done