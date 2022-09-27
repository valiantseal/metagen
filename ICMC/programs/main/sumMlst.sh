mkdir -p ./mlst_summary/

printf '%s,%s,%s\n' "sample" "organism" "ST" > ./mlst_summary/mlst_types.txt

for i in $(cat bacteria.list)
do

cd ./bactopia_gtdbtk/"$i"

bacteria=$(cat bacteria.id)

for sample in $(cat samples.list)
do
FILE=./output/"$sample"/mlst/default/ariba/mlst_report.tsv
if [ -f "$FILE" ]; then
    st=$(awk '{ print $1}' "$FILE" | sed -n '2p')
else 
    st='NONE'
fi

printf '%s,%s,%s\n' "$sample" "$bacteria" "$st" >> ../../mlst_summary/mlst_types.txt

done
cd ../../
done
