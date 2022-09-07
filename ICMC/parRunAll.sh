# version 1.01

#Need improvements
# write test for bactopia databasets
# write test to verify samples number for antimicrobial resistance
# rewrite checkGtdbTk.R so it would not be so messy and would include cases to deal with NA in metadata

sh ./programs/prepData.sh

cd process_par; ls -d */ | parallel -j 12 'cd {} && sh ../../programs/trimgalore.sh'; cd ../

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

Rscript --vanilla ./programs/findMismatchGtdk.R

bash -i ./programs/forBactopiaGtdbtk.sh

bash -i ./programs/bactopia.sh

echo "Finished Bactopia Run"
sleep 10s

bash -i ./programs/pangenome.sh

iqtree -s core_alignment.fasta -m TEST -nt 4 -b 200 # make a program to run with all samples

sh ./programs/pack.sh

sh ./programs/transferAws.sh









