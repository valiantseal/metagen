# c5.24xlarge
# 9.3h

par_samples=2

time sh ./programs/bin/prepSamples.sh # 1m45

cd gnuPar; time ls -d */ | parallel -j "$par_samples" 'cd {} && sh ../../programs/bin/filterMerge.sh'; cd ../ # 6m36

cd gnuPar; time ls -d */ | parallel -j "$par_samples" 'cd {} && sh ../../programs/bin/splitFiltered.sh'; cd ../ #1m16

time sh ./programs/bin/virema1.sh #485m22


time sh ./programs/bin/mergeVir1Out.sh #5m36

echo "virema1 merge is done"

cd gnuPar; time ls -d */ | parallel -j "$par_samples" 'cd {} && sh ../../programs/bin/prepVirema2.sh'; cd ../ # 7m13


time sh ./programs/bin/virema2Par.sh #75m36

