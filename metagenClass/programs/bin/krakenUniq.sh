conda activate krakenuniq

cd process

ls -d */ | parallel -j 4 'cd {} && krakenuniq --threads 4 --db /home/ubuntu/extraVol/krakenUniq/viral \
--report-file krakUniq_sample.report --output krakUniq_sample.kraken \
--classified-out krakUniq_classified.reads merged_reads.fa'



