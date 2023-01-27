conda activate /home/ubuntu/anaconda3/envs/seqtk

for i in $(cat newdir.list)
do
cd process/"$i"
seqtk subseq ./krakUniq_classified.reads \
krakenReads.id > selectKraken.reads
cd ../../
done
