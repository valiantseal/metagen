readsCh = Channel.fromFilePairs('/home/ubuntu/nextflowPipes/trial/input/*_R{1,2}.fastq.gz')

process trim {

  input: 
  tuple val(sample_id), path(sample_id)

  output:
  tuple val(sample_id), path('sample.fg.gz')

  script:
  """
  /home/ubuntu/trimGalore/TrimGalore-0.6.6/trim_galore --quality 20 \
  --paired $sample_id > ${sample_id}.fg.gz
  """
}

workflow {
    trim_ch = trim(readsCh)
    trim_ch.view()
}