for i in $(cat newdir.list); do mkdir -p ./trimKrUnVipr/"$i"/;
cp ./input/"$i"* ./trimKrUnVipr/"$i"/;
echo "$i" > ./trimKrUnVipr/"$i"/sample.name;
done