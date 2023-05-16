mkdir input


project="EHC_C19_2021"

dx select "$project"

dir="Ludy_Apr242023"

dx download "$dir"/*.fastq* -o ./input/
