library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

useCores <- 32

blastNames <- c('staxids', 'qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

filterReads <- function(df_path, blastNames) {
  x <- read_delim(df_path, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
  if (nrow(x) > 0) {
    colnames(x) <- blastNames
    x <- x[!is.na(x$bitscore), ]
    x <- x %>%
      filter(!grepl('Synthetic|clone', stitle, ignore.case = TRUE)) %>%
      arrange(qseqid, desc(bitscore)) %>%
      group_by(qseqid) %>%
      mutate(staxids = as.character(staxids)) %>%
      summarize(BID1 = first(staxids),
                BID2 = nth(unique(staxids), 2, default = NA_character_),
                BID3 = nth(unique(staxids), 3, default = NA_character_),
                .groups = 'drop') %>%
      ungroup()

    # Adjust for semicolons in taxonomic IDs
    x <- x %>%
      rowwise() %>%
      mutate(BID_combined = paste(BID1, BID2, BID3, sep = ";")) %>%
      mutate(BID_combined = str_replace_all(BID_combined, ";NA", ""),
             BID_combined = str_replace_all(BID_combined, "NA;", ""),
             BID_combined = str_replace_all(BID_combined, "NA", "")) %>%
      separate(BID_combined, into = c("BID1", "BID2", "BID3"), sep = ";", remove = FALSE, fill = "right") %>%
      select(-BID_combined)

    # Correct for duplicated or NA introductions in BID2 and BID3 due to shifting
    x <- x %>%
      mutate(BID2 = if_else(BID2 == BID1 | BID2 == "NA", NA_character_, BID2),
             BID3 = if_else(BID3 == BID1 | BID3 == BID2 | BID3 == "NA", NA_character_, BID3))

    return(x)
  } else {
    return(NULL)
  }
}

allDir <- list.files('./process', full.names = TRUE)

for (sampleDirPath in allDir) {
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)

  targDir <- paste0(sampleDirPath, '/splitSeq10K/')
  filesList <- list.files(targDir, full.names = TRUE)
  targSamples <- paste0(filesList, '/NtV4_blast.results')
  sampleName <- read.table(paste0(sampleDirPath, '/sample.name'), col.names = c("V1"))

  blastComb <- foreach(i = targSamples, .combine = rbind, .packages = c('readr', 'dplyr', 'stringr', 'tidyr')) %dopar% {
    filterReads(i, blastNames)
  }

  if (!is.null(blastComb)) {
    blastComb <- blastComb %>%
      rename(Read = qseqid) %>%
      select(Read, BID1, BID2, BID3) %>%
      mutate(Sample = sampleName$V1)

    outFile <- paste0(sampleDirPath, '/blastResTop_3.tsv')
    write.table(blastComb, file = outFile, col.names = TRUE, row.names = FALSE, quote = TRUE, sep = '\t')
  }

  stopCluster(cl)
  print(paste0(sampleDirPath, '________done!'))
}
