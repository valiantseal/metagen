setwd("/home/ubuntu/beast/Bri/50genes")

fasta<-phylotools::read.fasta('core50GenesAlign.fasta')
metaDat<-read.delim('first_ga-samples-with-date.txt',T, sep='\t')

fastaMet<-fasta[(fasta$seq.name%in%metaDat$sample),]

metaMiss<-metaDat[!(metaDat$sample%in%fastaMet$seq.name),]