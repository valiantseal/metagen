primers<-read.delim('./references/MPX_primers_220908.bed', F)

primers$V5<-0

primers$V6<-NA

for ( i in 1:nrow(primers)){
  if (grepl('_FOR', primers$V4[i])) {
    primers$V6[i]<-'+'
  } else if (grepl('_REV', primers$V4[i])) {
    primers$V6[i]<-'-'
  }
}

write.table(primers, './references/MPX_primers_220908_edit.bed', sep = '\t', row.names = F, col.names = F, quote = F)