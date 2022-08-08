/home/ubuntu/viral-ngs/databases/hg19.tar.gz
/home/ubuntu/viral-ngs/databases/metagenomics_contaminants_v3.tar.gz
/home/ubuntu/viral-ngs/databases/GRCh37.68_ncRNA-GRCh37.68_transcripts-HS_rRNA_mitRNA.tar.gz




bamDir=$(ls -d *_fastq_to_ubam)
bamFile=$(cat sample.txt).bam

miniwdl run \
https://raw.githubusercontent.com/broadinstitute/viral-ngs-staging/master/pipes/WDL/workflows/deplete_only.wdl \
deplete_taxa.raw_reads_unmapped_bam="$bamDir"/out/unmapped_bam/"$bamFile" \
deplete_taxa.bmtaggerDbs=/home/ubuntu/viral-ngs/databases/hg19.tar.gz \
deplete_taxa.bmtaggerDbs=/home/ubuntu/viral-ngs/databases/metagenomics_contaminants_v3.tar.gz \
deplete_taxa.bmtaggerDbs=/home/ubuntu/viral-ngs/databases/GRCh37.68_ncRNA-GRCh37.68_transcripts-HS_rRNA_mitRNA.tar.gz \
deplete_taxa.blastDbs=/home/ubuntu/viral-ngs/databases/metag_v3.ncRNA.mRNA.mitRNA.consensus.tar.gz




bamDir=$(ls -d *_fastq_to_ubam)
bamFile=$(cat sample.txt).bam
miniwdl run \
https://raw.githubusercontent.com/broadinstitute/viral-ngs-staging/master/pipes/WDL/workflows/assemble_denovo.wdl \
reads_unmapped_bam="$bamDir"/out/unmapped_bam/"$bamFile" \
reference_genome_fasta=/home/ubuntu/viral-ngs/MN908947.3Form.fasta \
trim_clip_db=/home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa \
deplete_bmtaggerDbs=/home/ubuntu/viral-ngs/databases/hg19.tar.gz \
deplete_bmtaggerDbs=/home/ubuntu/viral-ngs/databases/metagenomics_contaminants_v3.tar.gz \
deplete_bmtaggerDbs=/home/ubuntu/viral-ngs/databases/GRCh37.68_ncRNA-GRCh37.68_transcripts-HS_rRNA_mitRNA.tar.gz \
deplete_blastDbs=/home/ubuntu/viral-ngs/databases/metag_v3.ncRNA.mRNA.mitRNA.consensus.tar.gz \
call_isnvs=true


# fastp 
conda activate ivar
for i in $(cat sample.txt)
do 
r1=$(ls "$i"*R1*)
r2=$(ls "$i"*R2*)
fastp -i "$r1" -I "$r2" -o ./iSnvs/filt_"$r1" -O ./iSnvs/filt_"$r2" --detect_adapter_for_pe --thread 8
done

