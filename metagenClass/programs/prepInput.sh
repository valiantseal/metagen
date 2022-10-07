for i in $(cat newdir.list); do mkdir -p ./process/"$i"/;
cp ./input/"$i"* ./process/"$i"/;
echo "$i" > ./process/"$i"/sample.name;
done