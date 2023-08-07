
bamCompare -b1 ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-b2 ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam  \
-o figs/KJ-C1_scaleRc.bw \
--extendReads 150 \
--centerReads \
-p 12 \
--scaleFactorsMethod "readCount"

computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R macs2_poolAll/KJ-C_Chd7.bed \
-S figs/KJ-C1_scaleRc.bw \
--skipZeros \
-o figs/KJ-C1_scaleRc.mat \
-p 6 \
--outFileSortedRegions figs/regions_TSS_Chd7.bed

plotProfile -m figs/KJ-C1_scaleRc.mat \
-out figs/KJ-C1_scaleRc.png \
--refPointLabel "TSS" 

# sample and input separately
mkdir c1_figs

bamCoverage -b ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-o c1_figs/c1_bpmNorm.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 10

bamCoverage -b ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam \
-o c1_figs/i1_bpmNorm.bw \
--binSize 20 \
--normalizeUsing BPM \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 12

computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R macs2_poolAll/KJ-C_Chd7.bed \
-S c1_figs/*.bw \
--skipZeros \
-o c1_figs/c1_bp.mat \
-p 12 \
--outFileSortedRegions c1_figs/regions_TSS_Chd7.bed

plotProfile -m c1_figs/c1_bp.mat \
-out c1_figs/TSS_profile.png \
--perGroup \
--colors green purple \
--refPointLabel "TSS"

# try standard normalization
bamCoverage -b ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-o c1_figs/c1_rpkmNorm.bw \
--normalizeUsing "RPKM" \
--effectiveGenomeSize 2652783500 \
--extendReads 150 \
--centerReads \
-p 7

bamCoverage -b ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam \
-o c1_figs/i1_rpkmNorm.bw \
--normalizeUsing "RPKM" \
--effectiveGenomeSize 2652783500 \
--extendReads 150 \
--centerReads \
-p 7


computeMatrix reference-point --referencePoint TSS \
-b 1000 -a 1000 \
-R macs2_poolAll/KJ-C_Chd7.bed \
-S c1_figs/*_rpkmNorm.bw \
--skipZeros \
-o c1_figs/c1_rpkm.mat \
-p 12 \
--outFileSortedRegions c1_figs/regions_TSS_Chd7_peaks.bed

plotProfile -m c1_figs/c1_rpkm.mat \
-out c1_figs/TSS_profile_peaks_rpkm.png \
--perGroup \
--colors green purple \
--refPointLabel "TSS"

# try whole chd7 region
computeMatrix reference-point --referencePoint TSS \
-b 3000 -a 3000 \
-R Chd7.bed \
-S c1_figs/*_rpkmNorm.bw \
--skipZeros \
-o c1_figs/chd7_rpkm.mat \
-p 12 \
--outFileSortedRegions c1_figs/regions_TSS_Chd7_gene.bed

plotProfile -m c1_figs/chd7_rpkm.mat \
-out c1_figs/TSS_profile_chd7_rpkm.png \
--perGroup \
--colors green purple \
--refPointLabel "TSS"

# chd7 bpm
computeMatrix reference-point --referencePoint TSS \
-b 3000 -a 3000 \
-R Chd7.bed \
-S c1_figs/*_bpmNorm.bw \
--skipZeros \
-o c1_figs/chd7_bpm.mat \
-p 12 \
--outFileSortedRegions c1_figs/regions_TSS_Chd7_gene_bpm.bed

plotProfile -m c1_figs/chd7_bpm.mat \
-out c1_figs/TSS_profile_chd7_bpm.png \
--perGroup \
--colors green purple \
--refPointLabel "TSS"

# no normalization
bamCoverage -b ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-o c1_figs/c1_noNorm.bw \
--normalizeUsing "None" \
--extendReads 150 \
--centerReads \
-p 6

bamCoverage -b ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-o c1_figs/c1_noNorm.bw \
--normalizeUsing "None" \
--extendReads 150 \
--centerReads \
--binSize 1
-p 6

bamCoverage -b ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam \
-o c1_figs/i1_noNorm.bw \
--normalizeUsing "None" \
--extendReads 150 \
--centerReads \
-p 6

###
make_tracks_file --trackFiles macs2_poolAll/KJ-C_Chd7.bed c1_figs/c1_noNormBin1.bw  -o tracks.ini
pyGenomeTracks --tracks tracks.ini --region chr4:8848930-8854067 --outFileName nice_image.pdf

# process files with  RPKM and bin 1
bamCoverage -b ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam \
-o c1_figs/i1_rpkmNorm_bin1.bw \
--normalizeUsing "RPKM" \
--extendReads 150 \
--centerReads \
--binSize 1 \
--numberOfProcessors 6

bamCoverage -b ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-o c1_figs/c1_rpkmNorm_bin1.bw \
--normalizeUsing "RPKM" \
--extendReads 150 \
--centerReads \
--binSize 1 \
--numberOfProcessors 6


make_tracks_file --trackFiles macs2_poolAll/KJ-C_Chd7.bed c1_figs/c1_rpkmNorm_bin1.bw  -o tracks.ini
pyGenomeTracks --tracks tracks.ini --region chr4:8848930-8854067 --outFileName nice_image.png --height 35 --dpi 300

make_tracks_file --trackFiles macs2_poolAll/KJ-C_Chd7.bed c1_figs/i1_rpkmNorm_bin1.bw  -o tracks_i1.ini
pyGenomeTracks --tracks tracks_i1.ini --region chr4:8848930-8854067 --outFileName nice_image.png

make_tracks_file --trackFiles macs2_poolAll/KJ-C_Chd7.bed c1_figs/c1_rpkmNorm_bin1.bw macs2_poolAll/KJ-C_Chd7.bed c1_figs/i1_rpkmNorm_bin1.bw   -o tracks_comb.ini
pyGenomeTracks --tracks tracks_i1.ini --region chr4:8848930-8854067 --outFileName c1_i1_comb_pyGr.png --dpi 300
pyGenomeTracks --tracks tracks_comb.ini --region chr4:8848930-8854067 --outFileName c1_i1_comb_pyGr.png --dpi 300

# make c2 files
mkdir c2_figs
bamCoverage -b ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-o c2_figs/i2_rpkmNorm_bin1.bw \
--normalizeUsing "RPKM" \
--extendReads 150 \
--centerReads \
--binSize 1 \
--numberOfProcessors 6

bamCoverage -b ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-o c2_figs/c2_rpkmNorm_bin1.bw \
--normalizeUsing "RPKM" \
--extendReads 150 \
--centerReads \
--binSize 1 \
--numberOfProcessors 6

# adjust sampels by input
mkdir scale_input

bamCompare -b1 ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-b2 ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam  \
-o scale_input/c1_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

bamCompare -b1 ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-b2 ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-o scale_input/c2_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

# use merged files
bamCompare -b1 KJ-C.bam \
-b2 KJ-I.bam  \
-o scale_input/C_merged_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 12 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

# merged log2
bamCompare -b1 KJ-C.bam \
-b2 KJ-I.bam  \
-o scale_input/C_merged_scaleRc_log2.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 12 \
--scaleFactorsMethod "readCount" \
--operation "log2"

# H sampels
bamCompare -b1 ../output/bwa/mergedLibrary/KJ-H1_S17.mLb.clN.sorted.bam \
-b2 ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam  \
-o scale_input/h1_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

bamCompare -b1 ../output/bwa/mergedLibrary/KJ-H2_S18.mLb.clN.sorted.bam \
-b2 ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-o scale_input/h2_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"


# bowtie samples
# H
bamCompare -b1 ../output_bowtie2/bowtie2/mergedLibrary/KJ-H1_S17.mLb.clN.sorted.bam \
-b2 ../output_bowtie2/bowtie2/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam  \
-o bowtie2/h1_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

bamCompare -b1 ../output_bowtie2/bowtie2/mergedLibrary/KJ-H2_S18.mLb.clN.sorted.bam \
-b2 ../output_bowtie2/bowtie2/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-o bowtie2/h2_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

#C
bamCompare -b1 ../output_bowtie2/bowtie2/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-b2 ../output_bowtie2/bowtie2/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam  \
-o bowtie2/c1_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

bamCompare -b1 ../output_bowtie2/bowtie2/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-b2 ../output_bowtie2/bowtie2/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-o bowtie2/c2_scaleRc.bw \
--extendReads 150 \
--centerReads \
--numberOfProcessors 6 \
--scaleFactorsMethod "readCount" \
--operation "subtract"

