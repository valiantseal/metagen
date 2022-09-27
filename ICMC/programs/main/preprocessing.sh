aws s3 cp --recursive '9.15 - Witt_isolates' s3://transfer-files-emory/ICMC/Cecile/witt_isolates_09-15-22_fastq/


cd /home/ubuntu/ICMC/Cecile/witt_isolates_09-15-22
mkdir all_data
aws s3 cp --recursive s3://transfer-files-emory/ICMC/Cecile/witt_isolates_09-15-22_fastq/ ./all_data/

mkdir original_files
rm ./all_data/desktop.ini

mv ./all_data/*/*.fastq.gz ./original_files
