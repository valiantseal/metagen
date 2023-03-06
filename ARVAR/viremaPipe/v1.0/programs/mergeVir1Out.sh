
cd splitSeq

ls -d */ | parallel -j 35 'cd {} && samtools view -buSh -F 4 virema1.sam | samtools sort --threads 1 - -o virema1.bam'

cd ../


samtools merge finalBamFile.bam ./splitSeq/*/virema1.bam # extra copy to make sure that all files are merged

samtools merge virema1.bam ./splitSeq/*/virema1.bam