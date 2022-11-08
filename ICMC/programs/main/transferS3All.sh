# transfer
mkdir -p run_info
pwd > ./run_info/run.path
cat ./run_info/run.path |cut -f5 -d"/" > run_info/client.name
cat ./run_info/run.path |cut -f6 -d"/" > run_info/run.name

client=$(cat run_info/client.name)
run=$(cat run_info/run.name)

for i in $(cat bacteria.list)
do
aws s3 cp --recursive ./bactopia_gtdbtk/"$i"/output/ s3://transfer-files-emory/ICMC/"$client"/"$run"/"$i"/bactopia_output/
aws s3 cp --recursive ./bactopia_gtdbtk/"$i"/panOut/bactopia-tools/pangenome/pangenome/ s3://transfer-files-emory/ICMC/"$client"/"$run"/"$i"/pangenome/
done


echo s3://transfer-files-emory/ICMC/"$client"/"$run"/ > download.me
aws s3 cp download.me s3://transfer-files-emory/

