aws s3 cp --recursive "10.25 - K. pneumoniae FASTQ files/" s3://transfer-files-emory/ICMC/Ahmed/10.25-K.pneumoniae_fastq/


cd /home/ubuntu/ICMC/Ahmed
mkdir 10.25-K.pneumoniae

cd 
mkdir 10.25-K.pneumoniae

mkdir all_data

aws s3 cp --recursive s3://transfer-files-emory/ICMC/Ahmed/10.25-K.pneumoniae_fastq/ ./all_data/

mkdir original_files
rm ./all_data/desktop.ini

mv ./all_data/fastq/*.fastq.gz ./original_files
