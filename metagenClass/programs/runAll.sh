
sh ./programs/prepInput.sh
mkdir logs

bash -i ./programs/filterMerge.sh # fastp uses little bit more than 4 cpus max

bash -i ./programs/FqToFa.sh

bash -i ./programs/krakenUniq.sh

bash -i ./programs/splitReads.sh

sh ./programs/runBlastNt.sh

time sh ./programs/sumVirusBlastNt.sh

sh ./programs/getSelVirBlastNtRes.sh

Rscript --vanilla ./programs/countSummedReads.R

mkdir kraqSummary

for i in $(cat newdir.list);do cp process/"$i"/krakUniq_sample.report ./kraqSummary/"$i".report; done

