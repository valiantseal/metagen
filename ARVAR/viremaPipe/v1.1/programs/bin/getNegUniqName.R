#!/usr/bin/env Rscript


file1<-'samp_negative-mapped-reads_temp.txt'

file2<-'samp_positive-mapped-reads_temp.txt'

negative.list<-read.table(file1, F)

positive.list<-read.table(file2, F)

combined.list<-unique(rbind(negative.list, positive.list))

combined.sel<-combined.list[(combined.list$V1%in%negative.list$V1),]

outName<-"samp_negative-mapped-reads_unique_temp.txt"

write.table(combined.sel, outName, row.names = F, col.names = F, quote = F)