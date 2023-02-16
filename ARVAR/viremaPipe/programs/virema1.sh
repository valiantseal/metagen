cd splitSeq


ls -d */ | parallel -j 35 'cd {} && python /home/ubuntu/extraVol/Copyback/test_builds_2885y/programs/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna reads.fastq \
    virema1.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --p 1 -Overwrite >> ViReMaR1-results.txt; echo $(pwd) "____Virema_first_round_done"'
    