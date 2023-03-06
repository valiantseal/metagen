
for i in $(cat samples.list)
do
cd gnuPar/"$i"/splitSeq

ls -d */ | parallel -j 94 'cd {} && samtools view -buSh -F 4 virema1.sam | samtools sort --threads 1 - -o virema1.bam'

cd ../

samtools merge virema1.bam ./splitSeq/*/virema1.bam

cd  ../../

done
