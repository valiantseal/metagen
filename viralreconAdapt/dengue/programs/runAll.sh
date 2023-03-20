sh ./programs/download.sh

Rscript --vanilla ./programs/inputFile.R

sh ./programs/virDengue.sh

python ./programs/genomSum.py

python ./programs/results2S3.py
