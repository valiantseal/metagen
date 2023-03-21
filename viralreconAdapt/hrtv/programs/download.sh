mkdir input
mkdir water_input

project="HRTV"

dx select "$project"

#dir="2023-03-13"

fullPath=$(pwd)

dir=$(echo "$fullPath" | cut -f7 -d "/")

dx download "$dir"/1_raw_data/*.gz -o ./input/

mv ./input/Water* ./water_input/