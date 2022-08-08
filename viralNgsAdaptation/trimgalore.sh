conda activate ivar
for i in $(cat sample.txt)
do 
r1=$(ls "$i"*R1*)
r2=$(ls "$i"*R2*)
trim_galore --quality 20  --paired "$r1" "$r2" --nextera --output_dir ./iSnvs --cores 2
done