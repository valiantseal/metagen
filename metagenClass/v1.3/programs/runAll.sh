  GNU nano 6.2                                                                      runAll.sh                                                                               bash -i ./programs/bin/prepInput.sh

time bash -i ./programs/bin/filterMerge.sh

time bash -i ./programs/bin/FqToFa.sh

time bash -i ./programs/bin/krakenUniq.sh

cd process; time ls -d */ | parallel -j 4 'cd {} && Rscript --vanilla ../../programs/bin/sortKraken.R'; cd ..

time bash -i ./programs/bin/getSampKrakReads.sh

time bash -i ./programs/bin/splitReads.sh

time sh ./programs/bin/runBlastNt.sh

time Rscript --vanilla ./programs/bin/blastFiltTopSamp.R

time Rscript --vanilla ./programs/bin/confKrakBlast.R

nano output/term_mapping.txt

Rscript --vanilla ./programs/bin/confKrakBlast2.R
