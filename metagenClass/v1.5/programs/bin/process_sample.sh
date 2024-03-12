sample=$1
input_dir="$HOME/v1.5/input"
process_dir="$HOME/v1.5/process/$sample"
bbmap_dir="$HOME/bbmap"
trimmomatic_jar="$HOME/Trimmomatic-0.39/trimmomatic-0.39.jar"
adapters_path="$HOME/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa"
mkdir -p "$process_dir"
cp "$input_dir/${sample}_R1_001.fastq.gz" "$process_dir/"
cp "$input_dir/${sample}_R2_001.fastq.gz" "$process_dir/"
echo "$sample" > "$process_dir/sample.name"
R1="$process_dir/${sample}_R1_001.fastq.gz"
R2="$process_dir/${sample}_R2_001.fastq.gz"
R1_clumped="$process_dir/${sample}_R1_001_clumped.fastq.gz"
R2_clumped="$process_dir/${sample}_R2_001_clumped.fastq.gz"

"$bbmap_dir/clumpify.sh" in1="$R1" in2="$R2" out1="$R1_clumped" out2="$R2_clumped" dedupe

R1_paired="$process_dir/${sample}_R1_001_paired.fastq.gz"
R1_unpaired="$process_dir/${sample}_R1_001_unpaired.fastq.gz"
R2_paired="$process_dir/${sample}_R2_001_paired.fastq.gz"
R2_unpaired="$process_dir/${sample}_R2_001_unpaired.fastq.gz"

java -jar "$trimmomatic_jar" PE -phred33 "$R1_clumped" "$R2_clumped" "$R1_paired" "$R1_unpaired" "$R2_paired" "$R2_unpaired" ILLUMINACLIP:"$adapters_path":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Clean up
rm "$R1" "$R2" "$R1_clumped" "$R2_clumped" "$R1_unpaired" "$R2_unpaired"
mv "$R1_paired" "$process_dir/${sample}_R1_001.fastq.gz"
mv "$R2_paired" "$process_dir/${sample}_R2_001.fastq.gz"