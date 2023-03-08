for i in $(cat samples.list)
do
aws s3 cp --recursive ./bactopia_output/"$i"/ s3://transfer-files-emory/ICMC/Sarah/hflu_2023-02-27/original_bactopia_output/"$i"/
done

