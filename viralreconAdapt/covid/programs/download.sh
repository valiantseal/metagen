mkdir input
mkdir standard
mkdir noTrimOffset
mkdir -p water/input

conda activate dnaNexus
dx download 2022-09-05/1_raw_data/*fastq.gz -o ./input

mv ./input/Water* ./water/input/
