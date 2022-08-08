bamDir=$(ls -d *_align_and_plot)
bamFile=$(cat sample.txt).cleaned.mapped.bam

miniwdl run \
https://raw.githubusercontent.com/broadinstitute/viral-ngs-staging/master/pipes/WDL/workflows/isnvs_one_sample.wdl \
isnvs_per_sample.assembly_fasta=/home/ubuntu/viral-ngs/MN908947.3Form.fasta \
isnvs_per_sample.mapped_bam="$bamDir"//out/aligned_only_reads_bam/"$bamFile" 