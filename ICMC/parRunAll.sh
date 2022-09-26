# version 1.02
# run ec2 m5a.4xlarge

#Need improvements
# Fill the table for antimicrobial resistances with NONE if nothing found for the sample
# evaluate a test for bactopia databasets
# write test to verify samples number for pipeline outpusts: AMR, MLST, pangenome
# evaluate in whcih casses bactopia makes regular tree (not the fast tree) and edit so it would happen everytime
# for unresolved gtdb-tk taxa add editional tests with other methods kaiju/kraken/blast/gambit

sh ./programs/preprocessing.sh

sh ./programs/prepData.sh

cd process_par; ls -d */ | parallel -j 14 'cd {} && sh ../../programs/trimgalore.sh'; cd ../

echo "Finished Trimgalore"

# check number of created files output
ls process_par/*/*_val_1.fq.gz | wc -l; ls process_par/*/*_val_2.fq.gz | wc -l

sleep 10s


cd process_par; ls -d */ | parallel -j 2 'cd {} && sh ../../programs/parSpades.sh'; cd ../

echo "Finished Spades"
sleep 10s

sh ./programs/checkSpadesRes.sh

bash -i ./programs/gtdbtk.sh

echo "Finished GTDB-Tk"

Rscript --vanilla ./programs/formatMetadata.R

Rscript --vanilla ./programs/checkGtdb.R

bash -i ./programs/forBactopiaGtdbtk.sh

Rscript --vanilla ./programs/bactopDataTest.R

bash -i ./programs/bactopia.sh

echo "Finished Bactopia Run"
sleep 10s

bash -i ./programs/pangenome.sh

#iqtree -s core_alignment.fasta -m TEST -nt 4 -b 200 # may be make a program to run with all samples

sh ./programs/transferS3All.sh

# make summaries for cutsom output
sh ./programs/sumAMR.sh

bash -i ./programs/sumBactopia.sh

sh ./programs/sumMlst.sh

sh ./programs/transferS3Custom.sh












