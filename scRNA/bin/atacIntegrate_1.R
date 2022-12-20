library(Seurat)
library(Signac)
library(ggplot2)
library(GenomicRanges)
library(future)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

stress <- readRDS('atacAllFiltMacs2_stress')
control <- readRDS('atacAllFiltMacs2_control')

DefaultAssay(stress)<-'macs2'
DefaultAssay(control)<-'macs2'

stress_peak<-granges(stress)

control_peak<-granges(control)

identical(stress_peak, control_peak)

combined.peaks <- reduce(x = c(control_peak, stress_peak))

peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

control_frag<-Fragments(control)

stress_frag<-Fragments(stress)

identical(control_frag, stress_frag)

control_count<- FeatureMatrix(fragments = control_frag, features = combined.peaks)
#
stress_count<- FeatureMatrix(fragments = stress_frag, features = combined.peaks)

identical(control_count, stress_count)

###
control_chrom<-CreateChromatinAssay(control_count, fragments = control_frag)

stress_chrom<-CreateChromatinAssay(stress_count, fragments = stress_frag)

control[['Combined_peaks']]<-control_chrom
DefaultAssay(control)<-'Combined_peaks'

stress[['Combined_peaks']]<-stress_chrom
DefaultAssay(stress)<-'Combined_peaks'

control$dataset <- "Control"
stress$dataset <- "Stress"

saveRDS(stress, file = 'atacAllFiltMacs2_stress')

saveRDS(control, file = 'atacAllFiltMacs2_control')

rm(list = ls()[!ls() %in% c("control", "stress")])

gc()

