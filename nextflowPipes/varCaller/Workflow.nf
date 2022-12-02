nextflow.enable.dsl = 2
params.reads = './input/*_{1,2}.fastq.gz'
params.outdir = "./output"

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
    cpus 8

    /* 
       fastp process to remove adapters and low quality sequences
    */
    
    tag "filter $sample_id"
    
    publishDir "$params.outdir/${sample_id}/fastp/", 
    mode: 'copy'

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("filtered_{1,2}.fastq"), emit: filtered
    path("merged_prep_temp.fastq"), emit: merged
    path("${sample_id}.fastp.json"), emit: json
    path("virema1.fastq"), emit: virema1

 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w ${task.cpus} -j ${sample_id}.fastp.json \
      -m --merged_out merged_prep_temp.fastq -A -l 25 \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa;
      
      
    perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_1.fastq filtered_1_fastp-tagged_temp.fastq 1
    perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_2.fastq filtered_2_fastp-tagged_temp.fastq 2
    sh /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/editMerged.sh
    cat filtered_1_fastp-tagged_temp.fastq filtered_2_fastp-tagged_temp.fastq merged_prep_temp-tag.fastq >virema1.fastq
    
    """  
}  


workflow {
    fastp( reads )
}