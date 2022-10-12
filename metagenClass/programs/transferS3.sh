mkdir -p ./custom_output

cp blastKrakenConfirmedReadsTarget.csv ./custom_output/

cp blastNtSummary/blastNtSelVirSummary.csv ./custom_output/

cp kraqSummary/krakTargetVirSummary.csv ./custom_output/

aws s3 cp --recursive custom_output  s3://abombin/metagenClass/2022-02-25/custom_output/
