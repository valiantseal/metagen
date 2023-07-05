nextflow.enable.dsl = 2

params.reads = '/home/ubuntu/extraVol/hackaton2/data/fastqs/*_R{1,2}.fastq.gz'
params.kraken_db = "/home/ubuntu/extraVol/hackaton2/data"
params.amr_fasta = "/home/ubuntu/extraVol/hackaton2/data/AMR_CDS"
params.data_dir = "/home/ubuntu/extraVol/hackaton2/data"
params.outdir = "/home/ubuntu/extraVol/hackaton2/output"
params.amr_dir = "/home/ubuntu/extraVol/hackaton2/output/AMR_CDS"


// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


process index_amr {

    publishDir params.data_dir, mode:'copy'
    
    input:
    path fasta_file
    
    output:
    file "*"
    
    script:
    
    """
    kma index -i ${fasta_file} -o AMR_CDS
    
    """
  
}


process fastp {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filter $sample_id"

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}_filt_R{1,2}.fastq.gz")

 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_filt_R1.fastq.gz -O ${sample_id}_filt_R2.fastq.gz \
      --detect_adapter_for_pe -w 8 -j ${sample_id}.fastp.json \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa

    """  
}  

process run_kma {
    tag "AMR on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(filtered_reads)

    output:
    tuple val(sample_id), path("${sample_id}_kmamapped*")

    script:
    """
    
    kma -ipe  ${filtered_reads[0]}  ${filtered_reads[1]}  -t_db /home/ubuntu/extraVol/hackaton2/data/AMR_CDS -o ${sample_id}_kmamapped -mem_mode -ef -1t1 -cge -nf -vcf -t 1
    
    """
}


workflow {
  reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
  reference_fasta = Channel.fromPath(params.amr_fasta)
  index_amr( reference_fasta )
  filtered_ch = fastp( reads )
  kma_ch = run_kma(filtered_ch)
}