for i in $(cat samples.list)
do
cd gnuPar/"$i"

sh ~/extraVol/Copyback/test_builds_2885y/programs/virema2Par.sh

cd ../../

done
