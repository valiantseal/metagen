conda activate seqtk

for i in $(cat newdir.list)
do
cd process/"$i"
seqtk subseq ./krakUniq_classified.reads \
krakenReads.id > selectKraken.reads
cd ../../
done
