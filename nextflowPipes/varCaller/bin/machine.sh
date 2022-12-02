rm -f samples.list

#gunzip ./input/*.fastq.gz

for i in $(cd output && ls);
do
echo "$i" >> samples.list
machine=$(head -n 1 ./input/"$i"_1.fastq | sed 's/:/\t/g' - | awk '{print $1}' -)
echo "$machine" > ./output/"$i"/virema_1/machine.txt

done 
