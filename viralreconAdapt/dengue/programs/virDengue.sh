

/home/ubuntu/nextflow/nextflow run nf-core/viralrecon -r 2.5 \
    --max_cpus 30 \
    --max_memory '60.GB' \
    --input input.csv \
    --outdir output \
    --platform illumina \
    --protocol amplicon \
    --fasta ../../references/DENV1_reference_edit.fasta \
    --kraken2_db ../../references/kraken2-human-db \
    --skip_kraken2 \
    --save_reference false \
    --variant_caller 'ivar' \
    --primer_bed ../../references/DENV1_reference_edit.bed \
    --primer_left_suffix "_LEFT" \
    --primer_right_suffix "_RIGHT" \
    --ivar_trim_offset 5 \
    --skip_nextclade \
    --skip_assembly \
    -profile docker \
    -with-docker nfcore/virarecon


# run water

/home/ubuntu/nextflow/nextflow run nf-core/viralrecon -r 2.5 \
    --max_cpus 30 \
    --max_memory '60.GB' \
    --input inputWater.csv \
    --outdir water_output \
    --platform illumina \
    --protocol amplicon \
    --fasta ../../references/DENV1_reference_edit.fasta \
    --kraken2_db ../../references/kraken2-human-db \
    --skip_kraken2 \
    --save_reference false \
    --variant_caller 'ivar' \
    --primer_bed ../../references/DENV1_reference_edit.bed \
    --primer_left_suffix "_LEFT" \
    --primer_right_suffix "_RIGHT" \
    --ivar_trim_offset 5 \
    --skip_nextclade \
    --skip_assembly \
    -profile docker \
    -with-docker nfcore/virarecon

