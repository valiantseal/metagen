library(phylotools)

segments = c('HRTV_L', 'HRTV_M', 'HRTV_S')

for (segment in segments){
  inFasta = paste0('./references/', segment, '_consensus.fasta')
  fasta = phylotools::read.fasta(inFasta)
  fasta$seq.name = segment
  
  inBed = paste0('./references/', segment, '_primers.bed')
  primer = read.delim(inBed, F)
  primterFilt = primer[!(is.na(primer$V2)),]
  primterFilt$V4 = gsub(' ', '_', primterFilt$V4)
  primterFilt$V6[grepl('_REV', primterFilt$V1, ignore.case = T)]<-'-'
  primterFilt$V4[grepl('_REV', primterFilt$V1, ignore.case = T)]<-'Binding_region_REV'
  primterFilt$V6[grepl('_FWD', primterFilt$V1, ignore.case = T)]<-'+'
  primterFilt$V4[grepl('_FWD', primterFilt$V1, ignore.case = T)]<- 'Binding_region_FWD'
  primterFilt$V1 <- segment
  
  write.table(primterFilt, paste0('./references/', segment, '_primer_edit.bed'), sep = '\t', row.names = F, col.names = F, quote = F)
  phylotools::dat2fasta(fasta, paste0('./references/', segment, '_edit_consensus.fasta'))
}

