# m5a.8xlarge
/home/ubuntu/nextflow/nextflow run nf-core/viralrecon \
    --input input.csv \
    --outdir output \
    --platform illumina \
    --protocol amplicon \
    --genome 'MN908947.3' \
    --primer_bed /home/ubuntu/virilicon/AnneRegRuns/swift_refv3_primers.bed \
    --ivar_trim_offset 5 \
    --kraken2_db /home/ubuntu/kraken2/kraken2-human-db \
    --save_reference false \
    --variant_caller 'ivar' \
    --max_cpus 12 \
    --max_memory '50.GB' \
    --skip_assembly \
    -c /home/ubuntu/virilicon/AnneRegRuns/custom.config \
    -profile docker \
    -with-docker nfcore/virarecon | tee viralreconRun.log