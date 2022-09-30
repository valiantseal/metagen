run=${PWD##*/}

aws s3 cp --recursive output s3://transfer-files-emory/metagenClass/"$run"/


echo s3://transfer-files-emory/metagenClass/"$run"/ > download.me

aws s3 cp download.me s3://transfer-files-emory/
