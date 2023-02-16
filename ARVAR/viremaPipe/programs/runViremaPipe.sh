# took 20 hours to run

cd gnuPar; ls -d */ | parallel -j 2 'cd {} && sh ../../programs/viremaPipe.sh'; cd ../