for i in $(cat newdir.list)
do
cd process/"$i"/splitSeq10K

ls -d */ | parallel -j 28 'cd {} && sh ../../../../programs/bin/blastNtV4Par.sh'

cd ../../../
done