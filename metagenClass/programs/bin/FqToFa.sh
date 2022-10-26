conda activate seqtk

cd process

ls -d */ | parallel -j 4 'cd {} && seqtk seq -a merged_reads.fastq > merged_reads.fa'


