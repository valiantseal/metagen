# v1.2.1
# important fix added dependency program

#bash -i ./programs/bin/downloadDat.sh

bash -i ./programs/bin/downloadFromList.sh

bash -i ./programs/bin/prepInput.sh

time bash -i ./programs/bin/filterMerge.sh # changed to run in parallel

time bash -i ./programs/bin/FqToFa.sh 

time bash -i ./programs/bin/krakenUniq.sh

cd process; time ls -d */ | parallel -j 4 'cd {} && Rscript --vanilla ../../programs/bin/sortKraken.R'; cd ../ # 0.13

time bash -i ./programs/bin/getSampKrakReads.sh 

time bash -i ./programs/bin/splitReads.sh

time sh ./programs/bin/runBlastNt.sh # 26 min

#time Rscript --vanilla ~/extraVol/metagenClass/Dengue_virus/Batch_2/programs/bin/blastFiltCombSamp.R # 3.35, not a nessesary step for now

time Rscript --vanilla ./programs/bin/blastFiltTopSamp.R # 4.28

time Rscript --vanilla ./programs/bin/confKrakBlast.R # 0.14

sh ./programs/bin/resultToS3.sh

