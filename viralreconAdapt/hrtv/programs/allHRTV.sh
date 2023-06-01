
for i in $(cat programs/hrtv.segments)
do

/home/ubuntu/nextflow/nextflow run nf-core/viralrecon -r 2.5 \
    --max_cpus 30 \
    --max_memory '60.GB' \
    --input input.csv \
    --outdir output_"$i" \
    --platform illumina \
    --protocol amplicon \
    --fasta ../../references/"$i"_edit_consensus.fasta \
    --kraken2_db ../../references/kraken2-human-db \
    --skip_kraken2 \
    --save_reference false \
    --variant_caller 'ivar' \
    --primer_bed ../../references/"$i"_primer_edit.bed \
    --primer_left_suffix "_FWD" \
    --primer_right_suffix "_REV" \
    --ivar_trim_offset 5 \
    --skip_nextclade \
    --skip_assembly \
    -profile docker \
    -with-docker nfcore/virarecon

done

# run water
for i in $(cat programs/hrtv.segments)
do

/home/ubuntu/nextflow/nextflow run nf-core/viralrecon -r 2.5 \
    --max_cpus 30 \
    --max_memory '60.GB' \
    --input inputWater.csv \
    --outdir water_output_"$i" \
    --platform illumina \
    --protocol amplicon \
    --fasta ../../references/"$i"_edit_consensus.fasta \
    --kraken2_db ../../references/kraken2-human-db \
    --skip_kraken2 \
    --save_reference false \
    --variant_caller 'ivar' \
    --primer_bed ../../references/"$i"_primer_edit.bed \
    --primer_left_suffix "_FWD" \
    --primer_right_suffix "_REV" \
    --ivar_trim_offset 5 \
    --skip_nextclade \
    --skip_assembly \
    -profile docker \
    -with-docker nfcore/virarecon

done
