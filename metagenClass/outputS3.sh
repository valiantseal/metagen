run=${PWD##*/}

aws s3 cp --recursive output s3://transfer-files-emory/metagenClass/"$run"/

