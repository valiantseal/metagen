nextflow.enable.dsl = 2
params.reads = '/home/ubuntu/nextflow/trial/input/*_R{1,2}.fastq.gz'
params.outdir = "./denovo"

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


process fastp {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filter $sample_id"

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}_filt_R{1,2}.fastq.gz"), emit: filtered

 
    script:
    """
  echo sample_id
  /home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 \
  --paired ${reads[0]} ${reads[1]} --basename trimmed

  mv trimmed_val_1.fq.gz ${sample_id}_filt_R1.fastq.gz
  mv trimmed_val_2.fq.gz ${sample_id}_filt_R2.fastq.gz

    """  
}  

process assembly {
    /* 
       assembly step. here we used a conditional logic to choose the assembler
    */
    tag { sample_id }
    
    publishDir "$params.outdir/assemblies/", 
        mode: 'copy'
    
    input:
    tuple val(sample_id), path(filtered)  
    
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    /home/ubuntu/spades/SPAdes-3.15.4-Linux/bin/spades.py -1 ${filtered[0]} \
    -2 ${filtered[1]} -o spades_output \
    -t 8 --isolate
    mv ./spades_output/contigs.fasta ${sample_id}_contigs.fasta

    """
}

workflow {
    fastp( reads )
    assembly( fastp.out.filtered )
}