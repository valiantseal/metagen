
mkdir -p ./mlst_summary/

printf '%s,%s,%s\n' "sample" "organism" "ST" > ./mlst_summary/mlst_types.csv

for i in $(cat newdir.txt)
do

cd ./process/"$i"/shovill_out

#bacteria=$(cat bacteria.id)

FILE=./mlst.result
if [ -f "$FILE" ]; then
    st=$(awk '{ print $3}' "$FILE" | sed -n '1p')
    bacteria=$(awk '{ print $2}' "$FILE" | sed -n '1p')
else 
    st='NONE'
    bacteria='Not_Identified'
fi

printf '%s,%s,%s\n' "$i" "$bacteria" "$st" >> ../../../mlst_summary/mlst_types.csv

cd ../../../
done

