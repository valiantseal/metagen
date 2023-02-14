mkdir run_info
pwd > ./run_info/run.path
cat ./run_info/run.path |cut -f6 -d"/" > run_info/run.name

client="Viralrecon"
run=$(cat run_info/run.name)

mkdir output_files 
mv ./standard/output ./output_files/standard_output/
mv ./noTrimOffset/output_noTrimOffset/ ./output_files/output_noTrimOffset/
mv ./water/output/ ./output_files/output_water/

aws s3 cp --recursive ./output_files/ s3://transfer-files-emory/"$client"/"$run"/