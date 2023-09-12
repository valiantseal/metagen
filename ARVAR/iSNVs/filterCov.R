metaseq = read.csv("snvs_comb_res/metaseq_comb_derep_decont.csv")
ampseq = read.csv("snvs_comb_res/ampseq_comb_derep.csv")

cov_amp_new = read.csv("ampseqCovDepth.csv")
cov_amp_old = read.csv("ampseqOldSampCovDepth.csv")
comb_amp_cov = rbind(cov_amp_new, cov_amp_old)
colnames(comb_amp_cov)[1] = "OrigName"
comb_amp_cov = unique(comb_amp_cov)

cov_meta_new = read.csv("metaseqCovDepth.csv")
cov_meta_old = read.csv("metaseqOldSampCovDepth.csv")
comb_meta_cov = rbind(cov_meta_new, cov_meta_old)
colnames(comb_meta_cov)[1] = "OrigName"

metaseqCov = plyr::join(metaseq, comb_meta_cov, by="OrigName", type="left", match = "all")
any(is.na(metaseqCov$Coverage))

ampseqCov = plyr::join(ampseq, comb_amp_cov, by="OrigName", type="left", match = "all")
any(is.na(ampseqCov$Coverage))

metaseqCovFilt = metaseqCov[!metaseqCov$Coverage < 97,]
ampseqCovFit = ampseqCov[!ampseqCov$Coverage < 97,]

length(unique(metaseqCovFilt$Sample))
length(unique(ampseqCovFit$Sample))
length(unique(metaseqCovFilt$Sample[metaseqCovFilt$Sample%in%ampseqCovFit$Sample]))
overlapSamples = unique(metaseqCovFilt$Sample[metaseqCovFilt$Sample%in%ampseqCovFit$Sample])

write.csv(metaseqCovFilt, "snvs_comb_res/metaseq_comb_derep_decont_covFilt_97.csv", row.names = F)
write.csv(ampseqCovFit, "snvs_comb_res/ampseq_comb_derep_covFilt_97.csv", row.names = F)

# select only sample that overlap

metaOverlap = metaseqCovFilt[metaseqCovFilt$Sample%in%overlapSamples,]
ampOverlap = ampseqCovFit[ampseqCovFit$Sample%in%overlapSamples,]
length(unique(metaOverlap$Sample))
length(unique(ampOverlap$Sample))

write.csv(metaOverlap, "snvs_comb_res/metaseq_overlap_comb_derep_decont_covFilt_97.csv", row.names = F)
write.csv(ampOverlap, "snvs_comb_res/ampseq_overlap_comb_derep_covFilt_97.csv", row.names = F)

# get overlapping set for all coverage
overlapSamplesAll = unique(metaseqCov$Sample[metaseqCov$Sample%in%ampseqCov$Sample])

metaOverlapAll = metaseqCov[metaseqCov$Sample%in%overlapSamplesAll,]
ampOverlapAll = ampseqCov[ampseqCov$Sample%in%overlapSamplesAll,]

write.csv(metaOverlapAll, "snvs_comb_res/metaseq_overlap_comb_derep_decont_all.csv", row.names = F)
write.csv(ampOverlapAll, "snvs_comb_res/ampseq_overlap_comb_derep_all.csv", row.names = F)