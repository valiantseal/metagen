bamDir=$(ls -d *_deplete_only)
bamFile=$(cat sample.txt).cleaned.bam
assemblyDir=$(ls -d ../*_assemble_denovo)
assemblyFile=$(cat sample.txt).fasta

miniwdl run \
https://raw.githubusercontent.com/broadinstitute/viral-ngs-staging/master/pipes/WDL/workflows/align_and_plot.wdl \
align.reads_unmapped_bam="$bamDir"/out/cleaned_bam/"$bamFile" \
align.reference_fasta="$assemblyDir"/out/final_assembly_fasta/"$assemblyFile" \
align.aligner='bwa'

