mkdir -p gnuPar
rm -f samples.list
for i in $(cd input && ls *R1*) 
do
echo "$i" | cut -f1 -d "_" >> samples.list
done

for i in $(cat samples.list)
do
mkdir -p gnuPar/"$i"
cp ./input/"$i"* gnuPar/"$i"/
echo "$i" > gnuPar/"$i"/sample.name
done

cd gnuPar; ls -d */ | parallel -j 2 'cd {} && gzip -d *.gz'; cd ../