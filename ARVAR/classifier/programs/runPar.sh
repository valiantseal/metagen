
cd process; time ls -d */ | parallel -j 6 'cd {} && python ../../programs/drafts/bamReadCount.py --mapqual 0 --basequal 0 --allele_frequency 0 -r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../


python ../../programs/src/getVarFilesCLT.py -i {sample} --mapqual 5 --basequal 5 --allele_frequency 0 \
-r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 4

cd process_par; time ls -d */ | parallel -j 6 'cd {} && python ../../programs/drafts/bamReadCount.py --mapqual 0 --basequal 0 --allele_frequency 0 -r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../; sudo shutdown -h +10
cd amp_process_par; time ls -d */ | parallel -j 6 'cd {} && python ../../programs/drafts/bamReadCount.py \
--mapqual 0 --basequal 0 --allele_frequency 0 \
-r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 1'; cd ../; sudo shutdown -h +10

cd process_par; time ls -d */ | parallel -j 6 'cd {} && python ../../programs/drafts/lofreq.py \
--allele_frequency 0 -r ../../programs/data/MN908947.3.fna \
-b ../../programs/bin -t 1'; cd ../; sudo shutdown -h +10

