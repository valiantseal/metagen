setwd("/home/ubuntu/strain/ncov/data")
df<-read.table("MLS_all_meta.tsv", T, sep="\t")
df$date<-as.Date(df$date, format="%m/%d/%y")
df$date<-format(df$date, format="%Y-%m-%d")



df$date_submitted<-df$date
df$strain<-tolower(df$strain)

write.table(df, "MLS_all_metaForm.tsv", row.names = F, col.names = T, quote = F, sep = '\t')

fasta<-phylotools::read.fasta(file = "./MLS_all_seq.fasta", clean_name = F)

fasta$seq.name<-gsub("\\ .*","",fasta$seq.name)

phylotools::dat2fasta(fasta, outfile = "MLS_all_seqForm.fasta")