# v1.2.0

bash -i ./programs/bin/downloadDat.sh

bash -i ./programs/bin/prepInput.sh

time bash -i ./programs/bin/filterMerge.sh # changed to run in parallel

time bash -i ./programs/bin/FqToFa.sh 

time bash -i ./programs/bin/krakenUniq.sh


# wrewrting this part

cd process; time ls -d */ | parallel -j 4 'cd {} && Rscript --vanilla ../../programs/bin/sortKraken.R'; cd ../

###

time bash -i ./programs/bin/getSampKrakReads.sh 

time bash -i ./programs/bin/splitReads.sh

time sh ./programs/bin/runBlastNt.sh 

#time Rscript --vanilla ~/extraVol/metagenClass/Dengue_virus/Batch_2/programs/bin/blastFiltCombSamp.R # 3.35, not a nessesary step for now

time Rscript --vanilla ~/extraVol/metagenClass/Dengue_virus/Batch_2/programs/bin/blastFiltTopSamp.R # 4.28

#





time Rscript --vanilla ~/github/DailyWork/metagenClass/programs/development/combFilBlastAllSamp.R

#time Rscript --vanilla ./programs/bin/sumBlastReads.R

time Rscript --vanilla ./programs/bin/checkBlastKrakReadId.R


#cd testReads

#time sh ./programs/bin/transferS3.sh

time Rscript --vanilla ~/github/DailyWork/metagenClass/programs/development/filterVirReadsLength.R

time sh ~/github/DailyWork/metagenClass/programs/development/runGetFiltLenReads.sh # 40:51 processes, 57 min system time 4 times real time

time Rscript --vanilla ./programs/development/krakBlastTopMatch.R



