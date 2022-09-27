# version 1.03
# run ec2 m5a.4xlarge

#Need improvements
# Fill the table for antimicrobial resistances with NONE if nothing found for the sample
# evaluate a test for bactopia databasets
# write test to verify samples number for pipeline outpusts: AMR, MLST, pangenome
# evaluate in whcih casses bactopia makes regular tree (not the fast tree) and edit so it would happen everytime
# for unresolved gtdb-tk taxa add editional tests with other methods kaiju/kraken/blast/gambit

sh ./programs/main/preprocessing.sh

sh ./programs/main/prepData.sh

cd process_par; ls -d */ | parallel -j 7 'cd {} && sh ../../programs/main/trimgalore.sh'; cd ../

echo "Finished Trimgalore"

# check number of created files output
ls process_par/*/*_val_1.fq.gz | wc -l; ls process_par/*/*_val_2.fq.gz | wc -l

sleep 10s


cd process_par; ls -d */ | parallel -j 2 'cd {} && sh ../../programs/main/parSpades.sh'; cd ../

echo "Finished Spades"
sleep 10s

sh ./programs/tests/checkSpadesRes.sh

bash -i ./programs/main/gtdbtk.sh

echo "Finished GTDB-Tk"

Rscript --vanilla ./programs/main/formatMetadata.R

Rscript --vanilla ./programs/main/checkGtdb.R

bash -i ./programs/main/forBactopiaGtdbtk.sh

Rscript --vanilla ./programs/tests/bactopDataTest.R

bash -i ./programs/main/bactopia.sh

echo "Finished Bactopia Run"
sleep 10s

bash -i ./programs/main/pangenome.sh

#iqtree -s core_alignment.fasta -m TEST -nt 4 -b 200 # may be make a program to run with all samples

sh ./programs/main/transferS3All.sh

# make summaries for cutsom output
sh ./programs/main/sumAMR.sh

bash -i ./programs/main/sumBactopia.sh

sh ./programs/main/sumMlst.sh

sh ./programs/main/transferS3Custom.sh












