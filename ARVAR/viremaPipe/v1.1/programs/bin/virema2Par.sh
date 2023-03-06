for i in $(cat samples.list)
do
cd gnuPar/"$i"/


rm -rf splitSeq2

mkdir splitSeq2

lines=$(cat final_virema.fastq | wc -l)

parts=$(echo "$lines"/4/95+5|bc)
echo "$parts"

seqkit split final_virema.fastq -s "$parts" -j 4 -O ./splitSeq2/

cd splitSeq2

for j in *.fastq; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/final_virema.fastq; done

#ls -d */ | parallel -j 94 'cd {} && python /home/ubuntu/extraVol/Copyback/test_builds_2885y/viremaChuncks_1/programs/bin/ViReMa.py \
#    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna final_virema.fastq \
#    final_virema.sam \
#    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
#    --X 1 --MicroInDel_Length 2 --p 1 -Overwrite >> final_virema-results.txt; echo $(pwd) "____Virema_first_round_done"'
    

cd ../../../

    
done
