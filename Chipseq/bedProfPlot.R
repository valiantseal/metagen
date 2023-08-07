sample1 = read.delim("/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll/KJ-C.annotpeaks", T)
sampleSort = sample1[(sample1$Gene.Name == "Chd7"),]
sampleEdit = sampleSort[, 2:4]

write.table(sampleEdit, 'macs2_poolAll/KJ-C_Chd7.bed', row.names = F, col.names = F, quote = F, sep = "\t")