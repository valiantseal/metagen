params.index = "/home/ubuntu/nextflow/trial/input.csv"


process trim {

  input: 
  tuple val(sample), file(fastq_1), file(fastq_2)

  output:
  path("trimmed_val_{1,2}.fq.gz")

  script:
  """
  echo sample
  /home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 \
  --paired $fastq_1 $fastq_2 --basename trimmed
  """
}

process spades {
    input:
    path("x")

    output:
    path("spades_output/contigs.fasta")

    script:
    """
    echo $sample
    /home/ubuntu/spades/SPAdes-3.15.4-Linux/bin/spades.py -1 trimmed_val_1.fq.gz \
    -2 trimmed_val_2.fq.gz -o spades_output \
    -t 8 --isolate
    """
}

workflow {
    Channel.fromPath(params.index) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) } \
        | trim \
        | spades
}