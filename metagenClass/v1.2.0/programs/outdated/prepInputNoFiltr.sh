for i in $(cat newdir.list); do mkdir -p ./krUnVipr/"$i"/;
cp ./input/"$i"* ./krUnVipr/"$i"/;
echo "$i" > ./krUnVipr/"$i"/sample.name;
done