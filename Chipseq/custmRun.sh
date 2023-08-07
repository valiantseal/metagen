# run pooled controll and pooled all for cdh7 and check if my ranges overlap with the desired ones

time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C1_S15 -B --keep-dup all --extsize 200 --outdir macs2_poolContr

time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C2_S16 -B --keep-dup all --extsize 200 --outdir macs2_poolContr

mergePeaks -d given macs2_poolContr/KJ-C1_S15_peaks.narrowPeak macs2_poolContr/KJ-C2_S16_peaks.narrowPeak > macs2_poolContr/KJ-C_merge.bed
annotatePeaks.pl macs2_poolContr/KJ-C_merge.bed mm10  > macs2_poolContr/KJ-C_merge.annotpeaks

# running merge and annotation gives the same result on narrow peaks and bed files after ordering by start and end
mergePeaks -d given macs2_poolContr/KJ-C1_S15_peaks.bed macs2_poolContr/KJ-C2_S16_peaks.bed > macs2_poolContr/KJ-C_merge.bed
annotatePeaks.pl macs2_poolContr/KJ-C_merge.bed mm10  > macs2_poolContr/KJ-C_merge.annotpeaks

# all pooled peaks
mkdir macs2_poolAll

time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C -B --keep-dup all --extsize 200 --outdir macs2_poolAll

annotatePeaks.pl macs2_poolAll/KJ-C_peaks.narrowPeak mm10  > macs2_poolAll/KJ-C.annotpeaks

annotatePeaks.pl macs2_poolAll/KJ-C_peaks.bed mm10  > macs2_poolAll/KJ-C.annotpeaks

# pooled peaks with fdr 0.01
mkdir macs2_poolAll_fdr0.01
time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C -B --keep-dup all --extsize 200 --outdir macs2_poolAll_fdr0.01 -q 0.01

annotatePeaks.pl macs2_poolAll_fdr0.01/KJ-C_peaks.narrowPeak mm10  > macs2_poolAll_fdr0.01/KJ-C.annotpeaks


mkdir macs2_poolAll_fdr0.0001
time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C -B --keep-dup all --extsize 200 --outdir macs2_poolAll_fdr0.0001 -q 0.0001

annotatePeaks.pl macs2_poolAll_fdr0.0001/KJ-C_peaks.narrowPeak mm10  > macs2_poolAll_fdr0.0001/KJ-C.annotpeaks

# brod peaks
mkdir macs2_poolAll_fdr0.01_broad
time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C --keep-dup all --extsize 200 --broad --outdir macs2_poolAll_fdr0.01_broad -q 0.01

annotatePeaks.pl macs2_poolAll_fdr0.01_broad/KJ-C_peaks.broadPeak mm10  > macs2_poolAll_fdr0.01_broad/KJ-C.annotpeaks

# 0.05
mkdir macs2_poolAll_broad
time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-C --keep-dup all --extsize 200 --broad --outdir macs2_poolAll_broad -q 0.05

annotatePeaks.pl macs2_poolAll_broad/KJ-C_peaks.broadPeak mm10  > macs2_poolAll_broad/KJ-C.annotpeaks

# merge bam files
samtools merge -o KJ-C.bam ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam --threads 7
samtools merge -o KJ-I.bam ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam --threads 7
samtools index KJ-C.bam -@ 7
samtools index KJ-I.bam -@ 7

# make a run with homer to see the difference between homer and macs2
time makeTagDirectory Input_1_Tag/ ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam
time makeTagDirectory Input_2_Tag/ ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam

time makeTagDirectory C1_Tag/ ../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam
time makeTagDirectory C2_Tag/ ../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam

mkdir homer_merge
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -q 0.05 > homer_merge/homer_merged.annotpeaks

mkdir homer_merge_fold4
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 4 -q 0.05 > homer_merge_fold4/homer_merged.annotpeaks


mkdir homer_merge_fold1
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 1 -q 0.05 > homer_merge_fold1/homer_merged.annotpeaks

mkdir homer_merge_fold1_q0.1
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 1 -q 0.1 > homer_merge_fold1_q0.1/homer_merged.annotpeaks

# try to get them separately
mkdir homer_separate
findPeaks C1_Tag/ -style factor -F 2 -o homer_separate/C1.annotpeaks -i Input_1_Tag/
findPeaks C2_Tag/ -style factor -F 2 -o homer_separate/C2.annotpeaks -i Input_2_Tag/

annotatePeaks.pl homer_separate/C1.annotpeaks mm10  > homer_separate/C1_annot.annotpeaks
annotatePeaks.pl homer_separate/C2.annotpeaks mm10  > homer_separate/C2_annot.annotpeaks

mergePeaks homer_separate/C1.annotpeaks homer_separate/C2.annotpeaks > homer_separate/C_merged.peaks
annotatePeaks.pl homer_separate/C_merged.peaks mm10  > homer_separate/C_merged.annotpeaks


findPeaks C1_Tag/ -style histone -F 2 -o homer_separate/C1_hist.peaks -i Input_1_Tag/
findPeaks C2_Tag/ -style histone -F 2 -o homer_separate/C2_hist.peaks -i Input_2_Tag/

# macs2 narrow peaks hm
time macs2 callpeak -f AUTO -t ../output/bwa/mergedLibrary/KJ-H1_S17.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-H2_S18.mLb.clN.sorted.bam \
-c ../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam ../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam \
-g mm -n KJ-H -B --keep-dup all --extsize 200 --outdir macs2_poolAll_fdr0.01 -q 0.01

annotatePeaks.pl macs2_poolAll_fdr0.01/KJ-H_peaks.narrowPeak mm10  > macs2_poolAll_fdr0.01/KJ-H.annotpeaks

# hommer selected parameters
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 2 -q 0.05 -style histone > homer_merged_hist_fold2_fdr0.05.annotpeaks
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 2 -q 0.01 -style histone > homer_merged_hist_fold2_fdr0.01.annotpeaks
time getDifferentialPeaksReplicates.pl -t C1_Tag/ C2_Tag/ -i Input_1_Tag/ Input_2_Tag/ -genome mm10 -f 2 -q 0.01 -style factor > homer_merged_factor_fold2_fdr0.01.annotpeaks


# C histone
findPeaks C1_Tag/ -style histone -F 2 -fdr 0.01 -o homer_separate/C1_hist_fold2_fdr0.01.peaks -i Input_1_Tag/
findPeaks C2_Tag/ -style histone -F 2 -fdr 0.01  -o homer_separate/C2_hist_fold2_fdr0.01.peaks -i Input_2_Tag/
mergePeaks homer_separate/C1_hist_fold2_fdr0.01.peaks homer_separate/C2_hist_fold2_fdr0.01.peaks > homer_separate/C_hist_fold2_fdr0.01_merged.peaks
annotatePeaks.pl homer_separate/C_hist_fold2_fdr0.01_merged.peaks mm10  > homer_separate/C_hist_fold2_fdr0.01_merged.annotpeaks

findPeaks C1_Tag/ -style histone -F 4 -fdr 0.01 -o homer_separate/C1_hist_fold4_fdr0.01.peaks -i Input_1_Tag/
findPeaks C2_Tag/ -style histone -F 4 -fdr 0.01  -o homer_separate/C2_hist_fold4_fdr0.01.peaks -i Input_2_Tag/
mergePeaks homer_separate/C1_hist_fold4_fdr0.01.peaks homer_separate/C2_hist_fold4_fdr0.01.peaks > homer_separate/C_hist_fold4_fdr0.01_merged.peaks
annotatePeaks.pl homer_separate/C_hist_fold4_fdr0.01_merged.peaks mm10  > homer_separate/C_hist_fold4_fdr0.01_merged.annotpeaks

# C factor
findPeaks C1_Tag/ -style factor -F 2 -fdr 0.01 -o homer_separate/C1_factor_fold2_fdr0.01.peaks -i Input_1_Tag/
findPeaks C2_Tag/ -style factor -F 2 -fdr 0.01  -o homer_separate/C2_factor_fold2_fdr0.01.peaks -i Input_2_Tag/
mergePeaks homer_separate/C1_factor_fold2_fdr0.01.peaks homer_separate/C2_factor_fold2_fdr0.01.peaks > homer_separate/C_factor_fold2_fdr0.01_merged.peaks
annotatePeaks.pl homer_separate/C_factor_fold2_fdr0.01_merged.peaks mm10  > homer_separate/C_factor_fold2_fdr0.01_merged.annotpeaks

# reproduce prev results
findPeaks C1_Tag/ -style factor -F 2 -fdr 0.001 -o homer_separate/C1_factor_fold2_fdr0.001.peaks -i Input_1_Tag/
findPeaks C2_Tag/ -style factor -F 2 -fdr 0.001  -o homer_separate/C2_factor_fold2_fdr0.001.peaks -i Input_2_Tag/
mergePeaks homer_separate/C1_factor_fold2_fdr0.001.peaks homer_separate/C2_factor_fold2_fdr0.001.peaks > homer_separate/C_factor_fold2_fdr0.001_merged.peaks
annotatePeaks.pl homer_separate/C_factor_fold2_fdr0.001_merged.peaks mm10  > homer_separate/C_factor_fold2_fdr0.001_merged.annotpeaks

# homer find motifs
mkdir -p motifs/chd_all_shun; time findMotifsGenome.pl ../chd_all_shun.bed mm10 motifs/chd_all_shun -size 500 -mask
mkdir -p motifs/h3k_all_shun; time findMotifsGenome.pl ../h3k_all_shun.bed mm10 motifs/h3k_all_shun -size 500 -mask

# annotate peaks with mortifs homer
annotatePeaks.pl chd_chd7_shun.bed mm10 -m /home/flyhunter/miniconda3/envs/homer/share/homer/data/knownTFs/vertebrates/known.motifs > chd_chd7_shun_annotMotif.txt

