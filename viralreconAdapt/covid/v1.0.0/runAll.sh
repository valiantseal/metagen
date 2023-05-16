# write test to make sure that all input file samples have unique names

bash -i ./programs/download.sh

cd standard
Rscript --vanilla ../../programs/inputFile.R
sh ../../programs/viralrecon.sh

cd ../noTrimOffset
Rscript --vanilla ../../programs/inputFile.R
sh ../../programs/viralreconNoTrimOffset.sh

cd ../water
Rscript --vanilla ../../programs/inputFileWater.R
sh ../../programs/viralreconWater.sh

cd ../

sh ../programs/transferOutput.sh

