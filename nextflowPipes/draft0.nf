readsCh = Channel.fromFilePairs('/home/ubuntu/nextflow/trial/input/*_R{1,2}.fastq.gz')
params.outdir = 'my-results'

process trim {

  input: 
  tuple val(sample_id), path(sample_id)

  output:
  path "trimmed_*"

  script:
  """
  echo sample_id
  /home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 \
  --paired $sample_id --basename trimmed
  """
}

process spades {

    input:
    path x

    output:
    path "spades_output/contigs.fasta" 

    script:
    """
    /home/ubuntu/spades/SPAdes-3.15.4-Linux/bin/spades.py -1 trimmed_val_1.fq.gz \
    -2 trimmed_val_2.fq.gz -o spades_output \
    -t 8 --isolate
    """
}

process rename {
    publishDir "$params.outdir/contigs", mode: 'copy'
    input:
    path x

    output:
    path '*.fasta'

    script:
    """
    mv $x ${sample_id}_contigs.fasta
    """
}

workflow {
    trim_ch = trim(readsCh) | spades | rename
    trim_ch.view()
}
