nextflow.enable.dsl = 2
params.viremaSam = './output/*/virema_1/*_virema1.sam'
params.outdir = "./output"

virema1Sam = Channel.fromPath(params.viremaSam, checkIfExists: true).map {
    tuple( it.name.split('_')[0], it )
}
    
process samBam {

    cpus 8
    
    publishDir "$params.outdir/${sample_id}/virema_2/", 
    mode: 'copy'

    input:
    tuple val(sample_id), path(virema1Sam)

    output:
    tuple val(sample_id), path("${sample_id}_negative-strand.bam"), emit: negativeBam
    tuple val(sample_id), path("${sample_id}_positive-strand.bam"), emit: positiveBam
    tuple val(sample_id), path("${sample_id}_negative-mapped-reads_temp.txt"), emit: negativeMap
    tuple val(sample_id), path("${sample_id}_positive-mapped-reads_temp.txt"), emit: positiveMap
    tuple val(sample_id), path("${sample_id}_virema.fastq"), emit: forVirema2
    tuple val(sample_id), path("${sample_id}_virema-results.txt"), emit: viremaResult
    tuple val(sample_id), path("${sample_id}_virema.sam"), emit: viremaSam
    
    script:
    
    """
    samtools view -buSh -F 4 "${virema1Sam}" | samtools sort --threads ${task.cpus} - -o "${sample_id}_virema1.bam"

    samtools view -f 16 -b -o "${sample_id}_negative-strand.bam" "${sample_id}_virema1.bam"

    samtools view -F 16 -b -o "${sample_id}_positive-strand.bam" "${sample_id}_virema1.bam"
    
    samtools bam2fq "${sample_id}_negative-strand.bam" | grep -f /home/ubuntu/extraVol/Copyback/nextflowTrial/output/"${sample_id}"/virema_1/machine.txt - >"${sample_id}_negative-mapped-reads_temp.txt"
    
    samtools bam2fq "${sample_id}_positive-strand.bam" | grep -f /home/ubuntu/extraVol/Copyback/nextflowTrial/output/"${sample_id}"/virema_1/machine.txt - >"${sample_id}_positive-mapped-reads_temp.txt"
    
    Rscript --vanilla ~/extraVol/Copyback/nextflowTrial/bin/getNegUniqName.R
    
    sort -n "${sample_id}_positive-mapped-reads_temp.txt" | uniq -d - >> "${sample_id}_positive-mapped-reads_dedup_temp.txt"
    sort -n "${sample_id}_positive-mapped-reads_temp.txt" | uniq -u - >> "${sample_id}_positive-mapped-reads_dedup_temp.txt"
    
    sed "s/@//g" "${sample_id}_positive-mapped-reads_dedup_temp.txt" > "${sample_id}_positive-mapped-reads_dedup_temp1.txt"
    sed "s/@//g" "${sample_id}_negative-mapped-reads_unique_temp.txt" > "${sample_id}_negative-mapped-reads_unique_temp1.txt"

    seqtk subseq  /home/ubuntu/extraVol/Copyback/nextflowTrial/output/"${sample_id}"/fastp/virema1.fastq \
    "${sample_id}_negative-mapped-reads_unique_temp1.txt" > "${sample_id}_neg-mapped-reads_temp.fastq"


    seqtk subseq  /home/ubuntu/extraVol/Copyback/nextflowTrial/output/"${sample_id}"/fastp/virema1.fastq \
    "${sample_id}_positive-mapped-reads_dedup_temp1.txt" > "${sample_id}_pos-mapped-reads_temp.fastq"

    seqtk seq -r "${sample_id}_neg-mapped-reads_temp.fastq" > "${sample_id}_neg-mapped-rev_temp.fastq"

    cat "${sample_id}_neg-mapped-rev_temp.fastq" "${sample_id}_pos-mapped-reads_temp.fastq" > "${sample_id}_virema.fastq"
    
    python /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna "${sample_id}_virema.fastq" \
    "${sample_id}_virema.sam" \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p ${task.cpus} -Overwrite >>"${sample_id}_virema-results.txt"
    
    
    """
}

workflow {
    samBam( virema1Sam )
}
    