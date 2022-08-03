setwd("/home/ubuntu/beast/Bri/50genes")

df<-read.delim('/home/ubuntu/ICMC/Bri/bactopia_gtdbtk/haemophilus-influenzae/panOut/core_alignment_genePos.txt', F)

mergedIdSep <- data.frame(do.call('rbind', strsplit(as.character(df$V9),';',fixed=TRUE)))

data<-cbind(df, mergedIdSep)

dataFiltr<-data[!(data$X2=="gene=NA"),]

dataFiltr$length<-(dataFiltr$V5-dataFiltr$V4)

sum(dataFiltr$length)

#q1<-as.data.frame(dataFiltr[grepl("dacB", dataFiltr$X2),])

#q2<-as.data.frame(dataFiltr[grepl("16S", dataFiltr$X3),])

# select 100 random genes

set.seed(13) 
dfRand<-dataFiltr[sample(nrow(dataFiltr), 50), ]
sum(dfRand$length)



fastaSeq<-phylotools::read.fasta('/home/ubuntu/ICMC/Bri/bactopia_gtdbtk/haemophilus-influenzae/panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.fasta')

bySample<-split(fastaSeq, fastaSeq$seq.name)
sampleNames<-names(bySample)

startReg<-dfRand$V4
endReg<-dfRand$V5

sumData<-data.frame(matrix(ncol=2, nrow=0))
checkLen<-character()

for (i in sampleNames){
  df<-bySample[[i]]
  cutReg<-data.frame(substring(df$seq.text, first=startReg, last=endReg))
  colnames(cutReg)[1]<-'seq.txt'
  combReg<-data.frame(paste(cutReg$seq.txt, collapse = ''))
  checkLen<-append(checkLen, nchar(combReg[1,1]))
  sampleDf<-cbind(data.frame(i), combReg)
  colnames(sampleDf)<-names(fastaSeq)
  sumData<-rbind(sumData, sampleDf)
}

length(unique(checkLen))

phylotools::dat2fasta(sumData, outfile = "core50GenesAlign.fasta")