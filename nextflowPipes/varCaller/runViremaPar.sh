


rm -f ./gnuPar/*/*.txt
rm -f ./gnuPar/*/*.sam
rm -f ./gnuPar/*/*.bam
rm -f ./gnuPar/*/*.html
rm -f ./gnuPar/*/*.json
rm -f ./gnuPar/*/filtered*
rm -f ./gnuPar/*/merged*
rm -f ./gnuPar/*/merged*
rm -f ./gnuPar/*/virema*

cd gnuPar
ls -d */ | parallel -j 6 'cd {} && sh ~/extraVol/Copyback/nextflowTrial/viremaPipe.sh' 