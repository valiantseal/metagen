aws s3 cp --recursive 'EUHM CRE Outbreak Dec_Nov 2021' s3://transfer-files-emory/ICMC/Ahmed/euhm_cre_input/euhm_cre_2021/
aws s3 cp --recursive 'Order # 1001660 (2_3 rnds Rectal Swabs)' s3://transfer-files-emory/ICMC/Ahmed/euhm_cre_input/euhm_cre_2022/

cd /home/ubuntu/ICMC/Ahmed/euhm_cre
aws s3 cp --recursive s3://transfer-files-emory/ICMC/Ahmed/euhm_cre_input/ ./all_data/

mkdir original_files

mv '4313_S52_R2_001.fastq (1).gz' 4313_S52_R2_001.fastq.gz
mv '4314_S53_R1_001.fastq (1).gz' 4314_S53_R1_001.fastq.gz

mv ./all_data/euhm_cre_2021/*.fastq.gz ./original_files
mv ./all_data/euhm_cre_2022/*/*.fastq.gz ./original_files

aws s3 cp s3://transfer-files-emory/ICMC/Ahmed/euhm_cre_input/'Total_metadata_CRE_EUHM21_22_AB091922.xlsx' ./metadata/metadata.xlsx