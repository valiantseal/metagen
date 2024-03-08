library(tidyverse)
library(readr)

kraken_to_ncbi_mapping <- read_csv('../../kraken_to_ncbi_mapping.csv')
ncbi_parent_mapping <- setNames(kraken_to_ncbi_mapping$NCBI_Parent, kraken_to_ncbi_mapping$KrakenDB)


sampleName <- read.table('sample.name')
sampleName <- as.character(sampleName$V1)
virusId <- function() {

  df <- read.delim('krakUniq_sample.report', comment.char = "#")

  virus <- as.data.frame(df)
  colnames(virus)[colnames(virus) == "taxName"] <- "Virus" 
  virusID <- data.frame(unique(virus[, c('taxID', 'Virus')]))
  colnames(virusID)[1] <- 'krakTax' 
  return(virusID)
}

idList <- virusId()

readNameTab <- read.delim('krakUniq_sample.kraken', F, sep = '\t')
colnames(readNameTab)[3] <- 'krakTax'
colnames(readNameTab)[2] <- 'Read'
selCol <- readNameTab[, c('krakTax', 'Read')]
# join
cReads <- plyr::join(idList, selCol, by = 'krakTax', type = 'left', match = 'all')
allKrakenReads <- unique(cReads[, c('Virus', 'Read')])

allKrakenReads$Sample <- sampleName
write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

reads <- unique(allKrakenReads[, 'Read'])
write.table(reads, 'krakenReads.id', row.names = F, col.names = F, quote = F)

file_path <- "krakUniq_sample.kraken"

data <- read_lines(file_path)

parse_line <- function(line) {

    elements <- strsplit(line, "\\s+")[[1]]
    Read <- elements[2]
    primary_tax_id <- as.numeric(elements[3])

    tax_ids_all <- str_extract_all(line, "\\s[0-9]+(?=:|$)")[[1]]
    tax_ids_all <- as.numeric(gsub("\\s+", "", tax_ids_all))

    tax_ids_all <- sapply(tax_ids_all, function(id) {
        if (id > 1000000000) {
            parent_id <- ncbi_parent_mapping[as.character(id)]
            if (!is.na(parent_id) && parent_id <= 1000000000) {
                return(parent_id)
            }
        }
        return(id)
    })


    valid_tax_ids <- unique(tax_ids_all[tax_ids_all > 0])
    top_tax_ids <- head(valid_tax_ids, 3)

    if (length(top_tax_ids) < 3) {
        top_tax_ids <- c(top_tax_ids, rep(NA, 3 - length(top_tax_ids)))
    }


    return(data.frame(Read = Read, KID1 = top_tax_ids[1], KID2 = top_tax_ids[2], KID3 = top_tax_ids[3]))
}

remove_all_na_kids <- function(df) {
  df %>%
    filter(!is.na(KID1) | !is.na(KID2) | !is.na(KID3))
}

result <- bind_rows(lapply(data, parse_line)) %>%
  remove_all_na_kids()

write_csv(result, "krakenSelVirReads.tsv")
