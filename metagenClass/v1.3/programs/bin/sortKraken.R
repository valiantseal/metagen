# Read sample name
sampleName <- read.table('sample.name')
sampleName <- as.character(sampleName$V1)

virusId <- function() {
  # For each Kraken sample file, extract rows with potential viral sequences
  df <- read.delim('krakUniq_sample.report', comment.char = "#")
  # Make sure the "Virus" column contains taxonomic IDs (in this case "taxName")
  virus <- as.data.frame(df)
  colnames(virus)[colnames(virus) == "taxName"] <- "Virus"  # Rename the column
  virusID <- data.frame(unique(virus[, c('taxID', 'Virus')]))
  colnames(virusID)[1] <- 'krakId'
  return(virusID)
}

idList <- virusId()
# write.table(idList, 'kraken.ids', row.names = F, col.names = F, quote = F)
readNameTab <- read.delim('krakUniq_sample.kraken', F, sep = '\t')
colnames(readNameTab)[3] <- 'krakId'
colnames(readNameTab)[2] <- 'Read'
selCol <- readNameTab[, c('krakId', 'Read')]
# join
cReads <- plyr::join(idList, selCol, by = 'krakId', type = 'left', match = 'all')
allKrakenReads <- unique(cReads[, c('Virus', 'Read')])

# Note: Now "Virus" will contain the taxonomic IDs
allKrakenReads$Sample <- sampleName

# save sorted kraken output
write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

reads <- unique(allKrakenReads[, 'Read'])
write.table(reads, 'krakenReads.id', row.names = F, col.names = F, quote = F)
