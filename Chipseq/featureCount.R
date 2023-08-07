library(Rsubread)
library(Rsamtools)

df = read.delim("macs2_poolAll/KJ-C.annotpeaks")
colnames(df)[1] = "PeakID"

chd7 = df[(df$Gene.Name == "Chd7"),]

# check with macs2
#mcs2 = read.delim("macs2_poolAll/KJ-C_peaks.narrowPeak", F)
#macs2Sub = mcs2[(mcs2$V2 + 1) %in% chd7$Start,]

cur_names = c("GeneID", "Chr",	"Start",	"End",	"Strand")

chd_saf = chd7[, 1:5]
colnames(chd_saf) = cur_names

count_files = c("../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam", "../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam", 
                "../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam", "../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam")

chd7_counts = featureCounts(files = count_files, annot.ext = chd_saf, 
                            isGTFAnnotationFile = F, annot.inbuilt="mm10", nthreads = 7)

chd7_counts_df = data.frame(chd7_counts$counts)
chd7_counts_df$PeakID = rownames(chd7_counts_df)
colnames(chd7_counts_df)[1:4] = c("C1", "C2", "I1", "I2")


chd7_counts_df$C_sum = (chd7_counts_df$C1 + chd7_counts_df$C2) 
chd7_counts_df$I_sum = (chd7_counts_df$I1 + chd7_counts_df$I2)
sumFeatures = sum(chd7_counts_df$C_sum) + sum(chd7_counts_df$I_sum)

chd7_counts_df$C_peak_adj = chd7_counts_df$C_sum / sumFeatures
chd7_counts_df$I_peak_adj = chd7_counts_df$I_sum / sumFeatures
chd7_counts_df$Fold_change = chd7_counts_df$C_peak_adj  / chd7_counts_df$I_peak_adj
chd7_counts_df$Log2FCh = log2(chd7_counts_df$C_peak_adj  / chd7_counts_df$I_peak_adj)

shapiro.test(chd7_counts_df$Fold_change)
shapiro.test(chd7_counts_df$Log2FCh)

# adjust by the total number of reads in each library


c1_reads = countBam("../output/bwa/mergedLibrary/KJ-C1_S15.mLb.clN.sorted.bam")$records
c2_reads = countBam("../output/bwa/mergedLibrary/KJ-C2_S16.mLb.clN.sorted.bam")$records
c_reads_sum = c1_reads + c2_reads
i1_reads = countBam("../output/bwa/mergedLibrary/KJ-I1_S19.mLb.clN.sorted.bam")$records
i2_reads = countBam("../output/bwa/mergedLibrary/KJ-I2_S20.mLb.clN.sorted.bam")$records
i_reads_sum = i1_reads + i2_reads

chd7_counts_df$C_reads_adj = chd7_counts_df$C_sum / c_reads_sum
chd7_counts_df$I_reads_adj = chd7_counts_df$I_sum / i_reads_sum

chd7_counts_df$Fold_change_reads = chd7_counts_df$C_reads_adj / chd7_counts_df$I_reads_adj
chd7_counts_df$Log2FCh_reads = log2(chd7_counts_df$C_reads_adj / chd7_counts_df$I_reads_adj)

write.csv(chd7_counts_df, "macs2_poolAll/chd7_feature_counts.csv", row.names = F)

chd7CountSub = chd7_counts_df[, c("PeakID", "Fold_change", "Log2FCh", "Fold_change_reads", "Log2FCh_reads")]

combDat = plyr::join(chd7, chd7CountSub, by = "PeakID", type = "left", match = "all")


write.csv(combDat, "macs2_poolAll/chd7_peaks_annotatedfeature_counts.csv", row.names = F)


####
curDf = data.frame(chd7_counts$counts)
curDf $PeakID = rownames(curDf )
colnames(curDf )[1:4] = c("C1", "C2", "I1", "I2")

c1_sum = sum(curDf$C1)
c2_sum = sum(curDf$C2)
i1_sum = sum(curDf$I1)
i2_sum = sum(curDf$I2)

curDf$C1_adj = curDf$C1 / c1_sum 
curDf$C2_adj = curDf$C1 / c2_sum 
curDf$I1_adj = curDf$I1 / i1_sum
curDf$I2_adj = curDf$I2 / i2_sum 

curDf$C_mean =  (curDf$C1_adj + curDf$C2_adj) / 2
curDf$I_mean =  (curDf$I1_adj + curDf$I2_adj) / 2
curDf$Fold = curDf$C_mean / curDf$I_mean

curSub = curDf[, c("PeakID", "Fold")]

curCom = plyr::join(chd7, curSub, by = "PeakID", type = "left", match = "all")