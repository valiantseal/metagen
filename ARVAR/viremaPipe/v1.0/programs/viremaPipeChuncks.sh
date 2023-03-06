# c5.24xlarge
# 9.3h
cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh programs/filterMerge.sh'; cd ../

cd gnuPar; ls -d */ | parallel -j 2 'cd {} && sh programs/splitFiltered.sh'; cd ../

cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh programs/virema1.sh'; cd ../ #8.26h

cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh programs/mergeVir1Out.sh'; cd ../ 

echo "virema1 merge is done"

cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh programs/prepVirema2.sh'; cd ../ # 7.5m

#cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh ~/extraVol/Copyback/test_builds_2885y/programs/virema2.sh'; cd ../  # 3.75h

# alternatively to virema2.sh
# c5.9xlarge
cd gnuPar; time ls -d */ | parallel -j 2 'cd {} && sh programs/virema2Par.sh'; cd ../ # 78.19m

# check that working on two samples in parallel is faster than sequentially

#time sh programs/runVirema2Par.sh # 87 minutes

# taking into account such small difference between parallel and sequential samples for virema try to split it in more chuncks and run samples with for loop

