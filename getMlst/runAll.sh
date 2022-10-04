sh ./programs/prepInput.sh

cd process; ls -d */ | parallel -j 7 'cd {} && sh ../../programs/trimgalore.sh'; cd ../

# cd process; ls -d */ | parallel -j 2 'cd {} && bash -i ../../programs/shovill.sh'; cd ../ does not work

for i in $(cat newdir.txt); do cd process/"$i"/; bash -i ../../programs/shovill.sh; cd ../../; done

ls process/*/trimmed_val_1.fq.gz | wc -l; ls process/*/shovill_out/contigs.fa | wc -l

#mkdir -p mlst_summary

#printf '%s,%s,%s\n' "sample" "organism" "ST" > ./mlst_summary/mlst_types.tsv

bash -i ./programs/mlst.sh

sh ./programs/sumMlst.sh

#cd process; ls -d */ | parallel -j 13 'cd {} && bash -i ../../programs/mlst.sh'; cd ../