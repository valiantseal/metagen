# v1.00

bash -i ./programs/bin/downloadDat.sh

bash -i ./programs/bin/prepInput.sh

mkdir logs

time bash -i ./programs/bin/filterMerge.sh # 

time bash -i ./programs/bin/FqToFa.sh

time bash -i ./programs/bin/krakenUniq.sh

time bash -i ./programs/bin/splitReads.sh

time sh ./programs/bin/runBlastNt.sh # took more than 12 hours

time sh ./programs/bin/sumVirusBlastNt0.1.sh # took 101 minutes

time sh ./programs/bin/splitTargVirus.sh

Rscript --vanilla ./programs/bin/sumBlastReads.R

mkdir kraqSummary

for i in $(cat newdir.list);do cp process/"$i"/krakUniq_sample.report ./kraqSummary/"$i".report; done

Rscript --vanilla ./programs/bin/getKrakenIDs0.1.R

time sh ./programs/bin/getKrakSelReads.sh

Rscript --vanilla ./programs/bin/sumKrakReads.R

Rscript --vanilla ./programs/bin/checkBlastKrakReadId.R

sh ./programs/bin/transferS3.sh



#date > blast.time; time sh ./programs/bin/runBlastNt.sh; date >> blast.time; sudo shutdown -h +30




