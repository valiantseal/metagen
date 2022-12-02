nextflow.enable.dsl = 2
params.viremaSam = './output/*/virema_1/*_virema1.sam'
params.outdir = "./output"

virema1Sam = Channel.fromPath(params.viremaSam, checkIfExists: true).map {
    tuple( it.name.split('_')[0], it )
}
    
process samBam {

    cpus 4
    
    publishDir "$params.outdir/${sample_id}/samBam/", 
    mode: 'copy'

    input:
    tuple val(sample_id), path(virema1Sam)

    output:
    tuple val(sample_id), path("${sample_id}_negative-strand.bam"), emit: negativeBam
    tuple val(sample_id), path("${sample_id}_positive-strand.bam"), emit: positiveBam
    
    script:
    
    """
    samtools view -buSh -F 4 "${virema1Sam}" | samtools sort --threads ${task.cpus} - -o "${sample_id}_virema1.bam"

    samtools view -f 16 -b -o "${sample_id}_negative-strand.bam" "${sample_id}_virema1.bam"

    samtools view -F 16 -b -o "${sample_id}_positive-strand.bam" "${sample_id}_virema1.bam"
    
    """
}

workflow {
    samBam( virema1Sam )
}
    