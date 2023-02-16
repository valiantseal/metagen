

python /home/ubuntu/extraVol/Copyback/test_builds_2885y/programs/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna final_virema.fastq \
    final_virema.sam \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p 8 -Overwrite >> final_virema-results.txt