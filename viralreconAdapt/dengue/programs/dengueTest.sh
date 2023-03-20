
./nextflow run nf-core/viralrecon -r 2.5 \
    --max_cpus 16 \
    --max_memory '30.GB' \
    --input input.csv \
    --outdir output_dengue_test \
    --platform illumina \
    --protocol metagenomic \
    --fasta ./references/dengue2_ref.fasta \
    --kraken2_db ./references/kraken2-human-db \
    --save_reference false \
    --variant_caller 'bcftools' \
    --skip_assembly \
    --skip_nextclade \
    -profile docker \
    -with-docker nfcore/virarecon