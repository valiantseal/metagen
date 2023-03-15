for i in $(cat samples.list)
do
cd gnuPar/"$i"/splitSeq

ls -d */ | parallel -j 95 'cd {} && python ../../../../programs/bin/ViReMa.py \
    /home/ubuntu/references/MN908947.3.fna reads.fastq \
    virema1.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --p 1 -Overwrite >> ViReMaR1-results.txt; echo $(pwd) "____Virema_first_round_done"'

cd ../../../
done