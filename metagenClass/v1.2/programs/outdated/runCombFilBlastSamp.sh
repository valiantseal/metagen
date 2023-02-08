
cd ./process

ls -d */ | parallel -j 16 'cd {} && Rscript --vanilla ~/github/DailyWork/metagenClass/programs/development/combFilBlastSamp.R'

