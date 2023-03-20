mkdir input
mkdir water_input

project="Dengue_virus"

dx select "$project"

dir="2023-03-14"

dx download "$dir"/1_raw_data/*.gz -o ./input/

mv ./input/Water* ./water_input/