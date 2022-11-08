rm ./kraqSummary/*.readnames
rm ./kraqSummary/*.reads

cd process
ls -d */ | parallel -j 16 'cd {} && grep -F -f  ../../kraqSummary/kraken.ids ./krakUniq_sample.kraken > selKraken.readnames'
cd ../

for sample in $(cat newdir.list);
do
cat ./process/"$sample"/selKraken.readnames >> ./kraqSummary/"$sample".reads
done