conda activate krakenuniq

cd process

ls -d */ | parallel -j 5 'cd {} && krakenuniq --threads 5 --db /home/ubuntu/extraVol/krakenUniq/viral \
--report-file krakUniq_sample.report --output krakUniq_sample.kraken \
--classified-out krakUniq_classified.reads merged_reads.fa'



