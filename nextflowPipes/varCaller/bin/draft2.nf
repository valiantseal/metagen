nextflow.enable.dsl = 2
params.reads = './input/*_R{1,2}.fastq'
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
    tuple val(sample_id), path("${sample_id}_virema-results.txt"), emit: viremaResult
    tuple val(sample_id), path("${sample_id}_virema.sam"), emit: viremaSam
    

 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o filtered_1.fastq -O filtered_2.fastq \
      --detect_adapter_for_pe -w ${task.cpus} -j ${sample_id}.fastp.json \
      -m --merged_out merged_prep_temp.fastq -A -l 25 \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa;
      
      
    perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_1.fastq filtered_1_fastp-tagged_temp.fastq 1 & perl /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/AddPairedEndSuffix.pl filtered_2.fastq filtered_2_fastp-tagged_temp.fastq 2;
    sh /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/editMerged.sh
    cat filtered_1_fastp-tagged_temp.fastq filtered_2_fastp-tagged_temp.fastq merged_prep_temp-tag.fastq >virema1.fastq
    
    
    python /home/ubuntu/extraVol/Copyback/nextflowTrial/bin/ViReMa.py \
    /home/ubuntu/extraVol/Copyback/nextflowTrial/reference/GCA_009858895.3_ASM985889v3_genomic.200409.fna virema1.fastq \
    "${sample_id}_virema1.sam" \
    --Output_Dir ./ --Seed 13 --ErrorDensity 2,20 \
    --X 1 --MicroInDel_Length 2 --Chunk 10000000 --p ${task.cpus} -Overwrite >>"${sample_id}_ViReMaR1-results.txt"


    samtools view -buSh -F 4 "${sample_id}_virema1.sam" | samtools sort --threads ${task.cpus} - -o "${sample_id}_virema1.bam"
    samtools view -f 16 -b -o "${sample_id}_negative-strand.bam" "${sample_id}_virema1.bam"
    samtools view -F 16 -b -o "${sample_id}_positive-strand.bam" "${sample_id}_virema1.bam"
    
    samtools bam2fq "${sample_id}_negative-strand.bam" | grep -f /home/ubuntu/extraVol/Copyback/nextflowTrial/machine_dir/"${sample_id}_machine.txt" - >"${sample_id}_negative-mapped-reads_temp.txt"
    samtools bam2fq "${sample_id}_positive-strand.bam" | grep -f /home/ubuntu/extraVol/Copyback/nextflowTrial/machine_dir/"${sample_id}_machine.txt" - >"${sample_id}_positive-mapped-reads_temp.txt"
    Rscript --vanilla ~/extraVol/Copyback/nextflowTrial/bin/getNegUniqName.R
    
    sort -n "${sample_id}_positive-mapped-reads_temp.txt" | uniq -d - >> "${sample_id}_positive-mapped-reads_dedup_temp.txt"
    sort -n "${sample_id}_positive-mapped-reads_temp.txt" | uniq -u - >> "${sample_id}_positive-mapped-reads_dedup_temp.txt"
    sed "s/@//g" "${sample_id}_positive-mapped-reads_dedup_temp.txt" > "${sample_id}_positive-mapped-reads_dedup_temp1.txt"
    sed "s/@//g" "${sample_id}_negative-mapped-reads_unique_temp.txt" > "${sample_id}_negative-mapped-reads_unique_temp1.txt"
    
    seqtk subseq virema1.fastq "${sample_id}_negative-mapped-reads_unique_temp1.txt" > "${sample_id}_neg-mapped-reads_temp.fastq"
    seqtk subseq virema1.fastq "${sample_id}_positive-mapped-reads_dedup_temp1.txt" > "${sample_id}_pos-mapped-reads_temp.fastq"
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
    fastp( reads )
    
}