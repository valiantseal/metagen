conda activate ivar
for i in $(cat sample.txt)
do 
r1=$(ls "$i"*R1*)
r2=$(ls "$i"*R2*)
fastp -i "$r1" -I "$r2" -o ./iSnvs/filt_"$r1" -O ./iSnvs/filt_"$r2" --detect_adapter_for_pe \
--adapter_sequence /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa \
--adapter_sequence_r2 /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa --thread 8
done