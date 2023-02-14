mkdir input
mkdir water_input

project="Mpox"

dir=$(basename "$PWD") 

dx select "$project"

dx download "$dir"/1_raw_data/*.fastq.gz -o ./input

mv ./input/Water* ./water_input/