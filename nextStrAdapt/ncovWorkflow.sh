cd /home/ubuntu/strain/processData/georgia/2021-03-01_2022-06-23
aws s3 cp s3://transfer-files-emory/virus/Anne/data/georgia/gisaidAuspiceInputHcov19Georgia_2021-03-01_2022-06-23.zip ./
for f in *.tar; do tar xf "$f"; done

ls *.tar | wc -l
ls *.fasta | wc -l
ls *.tsv | wc -l

rm *.tar

mkdir joined
Rscrtipt -- vanilla ../../../../programs/joinData.R
cp ./joined/* /home/ubuntu/strain/ncov/data/

cd /home/ubuntu/strain/ncov/my_profiles
mkdir georgia20220623Global
cd /home/ubuntu/strain/ncov
snakemake --profile my_profiles/georgia20220623Global -p