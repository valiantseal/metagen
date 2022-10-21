readsCh = Channel.fromFilePairs('/home/ubuntu/nextflow/trial/input/*_R{1,2}.fastq.gz')

process trim {

  input: 
  tuple val(sample_id), path(sample_id)

  output:
  tuple val(sample_id), path("trimmed_val_{1,2}.fq.gz")

  script:
  """
  echo sample_id
  /home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 \
  --paired $sample_id --basename trimmed
  """
}

workflow {
    trim_ch = trim(readsCh)
    trim_ch.view()
}
