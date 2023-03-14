for i in ls -d *output*
do
aws s3 cp --recursive "$i"/ s3://transfer-files-emory/Viralrecon/HRTV/2023-02-23/"$i"/
done
