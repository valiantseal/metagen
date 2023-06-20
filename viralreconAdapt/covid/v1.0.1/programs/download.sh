mkdir input


project="EHC_C19_2021"

dx select "$project"

dir="for_andrei_june152023/1_raw_data"

dx download "$dir"/*.fastq* -o ./input/
