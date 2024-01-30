library(tidyverse)
library(readr)

kraken_to_ncbi_mapping <- read_csv('../../kraken_to_ncbi_mapping.csv')
ncbi_parent_mapping <- setNames(kraken_to_ncbi_mapping$NCBI_Parent, kraken_to_ncbi_mapping$KrakenDB)

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
  colnames(virusID)[1] <- 'krakTax'  # Changed from 'krakId' to 'krakTax'
  return(virusID)
}

idList <- virusId()
# write.table(idList, 'kraken.ids', row.names = F, col.names = F, quote = F)
readNameTab <- read.delim('krakUniq_sample.kraken', F, sep = '\t')
colnames(readNameTab)[3] <- 'krakTax'  # Changed from 'krakId' to 'krakTax'
colnames(readNameTab)[2] <- 'Read'  # Capitalize 'read' to 'Read'
selCol <- readNameTab[, c('krakTax', 'Read')]
# join
cReads <- plyr::join(idList, selCol, by = 'krakTax', type = 'left', match = 'all')  # Changed from 'krakId' to 'krakTax'
allKrakenReads <- unique(cReads[, c('Virus', 'Read')])

# Note: Now "Virus" will contain the taxonomic IDs
allKrakenReads$Sample <- sampleName

# save sorted kraken output
write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

reads <- unique(allKrakenReads[, 'Read'])
write.table(reads, 'krakenReads.id', row.names = F, col.names = F, quote = F)

# Define the path to the file
file_path <- "krakUniq_sample.kraken"

# Read the file line by line
data <- read_lines(file_path)

parse_line <- function(line) {
    # Split the line into elements
    elements <- strsplit(line, "\\s+")[[1]]
    # Extract the Read and the primary taxonomic ID
    Read <- elements[2]
    primary_tax_id <- as.numeric(elements[3])

    # Extract all potential tax IDs from the line
    tax_ids_all <- str_extract_all(line, "\\s[0-9]+(?=:|$)")[[1]]
    tax_ids_all <- as.numeric(gsub("\\s+", "", tax_ids_all))

# Replace tax IDs over 1,000,000,000 with their NCBI parents
    tax_ids_all <- sapply(tax_ids_all, function(id) {
        if (id > 1000000000) {
            parent_id <- ncbi_parent_mapping[as.character(id)]
            if (!is.na(parent_id) && parent_id <= 1000000000) {
                return(parent_id)
            }
        }
        return(id)
    })

    # Take the first three unique valid tax IDs
    valid_tax_ids <- unique(tax_ids_all[tax_ids_all > 0])
    top_tax_ids <- head(valid_tax_ids, 3)

    # Pad the result with NA if there are less than 3 tax IDs
    if (length(top_tax_ids) < 3) {
        top_tax_ids <- c(top_tax_ids, rep(NA, 3 - length(top_tax_ids)))
    }

    # Return the Read and the valid taxonomic IDs
    return(data.frame(Read = Read, KID1 = top_tax_ids[1], KID2 = top_tax_ids[2], KID3 = top_tax_ids[3]))
}

# Apply the function to each line and combine results into a dataframe
result <- bind_rows(lapply(data, parse_line))

# Write to CSV
write_csv(result, "krakenSelVirReads.tsv")
