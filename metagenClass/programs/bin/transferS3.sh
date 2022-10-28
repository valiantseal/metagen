mkdir -p ./custom_output

cp blastKrakenConfirmedReadsTarget.csv ./custom_output/

cp blastNtSummary/blastNtSelVirSummary.csv ./custom_output/

cp kraqSummary/krakTargetVirSummary.csv ./custom_output/

aws s3 cp --recursive custom_output  s3://abombin/metagenClass/Nima_p2/custom_output/


cd testReads
for i in $(ls *.fasta); do aws s3 cp "$i" s3://abombin/metagenClass/Nima_p2/custom_output/confirmedReads/"$i"; done
