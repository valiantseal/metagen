# v1.01

bash -i ./programs/bin/downloadDat.sh

bash -i ./programs/bin/prepInput.sh

mkdir logs

time bash -i ./programs/bin/filterMerge.sh # really can use 16 cores per sample, dows not want to work with parallel

time bash -i ./programs/bin/FqToFa.sh 

time bash -i ./programs/bin/krakenUniq.sh

for i in $(cat newdir.list);do cp process/"$i"/krakUniq_sample.report ./kraqSummary/"$i".report; done

Rscript --vanilla ./programs/bin/getKrakenIDs0.1.R

time sh ./programs/bin/getKrakSelReads0.1.sh

time Rscript --vanilla ./programs/bin/sumKrakReads.R

time Rscript --vanilla ./programs/bin/localKrakReads.R

time bash -i ./programs/bin/getSampKrakReads.sh

time bash -i ./programs/bin/splitReads.sh

time sh ./programs/bin/runBlastNt.sh 

time sh ./programs/bin/sumVirusBlastNt0.1.sh 

time sh ./programs/bin/splitTargVirus.sh

time Rscript --vanilla ./programs/bin/sumBlastReads.R

time Rscript --vanilla ./programs/bin/checkBlastKrakReadId.R

mkdir testReads

cd testReads


time sh ./programs/bin/transferS3.sh



#date > blast.time; time sh ./programs/bin/runBlastNt.sh; date >> blast.time; sudo shutdown -h +30



