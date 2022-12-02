setwd("/home/ubuntu/strain/ncov/data")

mlsMeta<-read.delim("MLS_teammates_meta.tsv", T, sep = '\t')

mlsFasta<-phylotools::read.fasta("MLS_teammates_seq.fasta")

namesCheck<-mlsFasta[(mlsFasta$seq.name%in%mlsMeta$strain),]

refMeta<-read.delim('USA_AY4_meta.tsv', T, sep='\t')

mlsMeta$location[mlsMeta$location==""]<-"NotMls"

refMeta$location[refMeta$location==""]<-"NotMls"

write.table(mlsMeta, "MLS_teammates_metaEdit.tsv", row.names = F, quote = F, sep = '\t')
write.table(refMeta, "USA_AY4_metaEdit.tsv", row.names = F, quote = F, sep = '\t')


## edit for the AY4all

ayMeta<-read.delim("AY4all.meta.tsv", T, sep='\t')

ayMeta$date<-as.Date(ayMeta$date, format="%m/%d/%y")
ayMeta$date<-format(ayMeta$date, format="%Y-%m-%d")

ayFasta<-phylotools::read.fasta('AY4all.sequences.fasta')


ayMetaSel<-ayMeta[(ayMeta$strain%in%ayFasta$seq.name),]

write.table(ayMeta, "AY4all_metaEdit.tsv", row.names = F, col.names = T, quote = F, sep = '\t')
