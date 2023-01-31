conda activate /home/ubuntu/anaconda3/envs/dnaNexus

project="BWH_mNGS"

dir=$(basename "$PWD") 


mkdir input
mkdir process

dx select "$project"


#for i in $(cat downloadSamples.list)
#do
#dx download "$i"* -o ./input
#done

dx download "$dir"/1_raw_data/*.fastq.gz -o ./input