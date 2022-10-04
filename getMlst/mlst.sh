conda activate shovill

#mkdir -p mlst_summary

#printf '%s,%s,%s\n' "sample" "organism" "ST" > ./mlst_summary/mlst_types.tsv


for i in $(cat newdir.txt)
do
cd process/"$i"/shovill_out
mlst contigs.fa > mlst.result
#cut -f1,2,3 mlst.result >> ../../../mlst_summary/mlst_types.tsv
cd ../../../
done

#sed 's/\t/,/g' ./mlst_summary/mlst_types.tsv > ./mlst_summary/mlst_types.csv

#cd shovill_out
#mlst contigs.fa > mlst.result
#cut -f1,2,3 mlst.result >> ../../../mlst_summary/mlst_types.tsv
