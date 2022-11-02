conda activate /home/ubuntu/anaconda3/envs/seqtk

cd process

ls -d */ | parallel -j 16 'cd {} && seqtk seq -a merged_reads.fastq > merged_reads.fa'


