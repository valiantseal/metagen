curDir=$(pwd)

project=$(echo "$curDir" |cut -f5 -d "/")

division=$(echo "$curDir" |cut -f6 -d "/")

run=$(echo "$curDir" |cut -f7 -d "/")

aws s3 cp --recursive ./main_output s3://transfer-files-emory/"$project"/"$division"/"$run"/main_output/

aws s3 cp --recursive ./water_output s3://transfer-files-emory/"$project"/"$division"/"$run"/water_output/