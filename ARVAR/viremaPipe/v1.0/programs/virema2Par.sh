rm -rf splitSeq2

mkdir splitSeq2

seqkit split final_virema.fastq -s 1000000 -j 4 -O ./splitSeq2/

cd splitSeq2

for j in *.fastq; do mkdir "$j"_dir; mv "$j" ./"$j"_dir/final_virema.fastq; done

ls -d */ | parallel -j 16 'cd {} && python /home/ubuntu/extraVol/Copyback/test_builds_2885y/programs/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna final_virema.fastq \
    final_virema.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p 1 -Overwrite >> final_virema-results.txt; echo $(pwd) "____Virema_first_round_done"'
    

