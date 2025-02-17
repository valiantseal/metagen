library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(stringr)

useCores <- 15

blastNames <- c('staxids', 'qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

filterReads <- function(df) {
  x <- readr::read_delim(df, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
  if (nrow(x) > 0) {
    colnames(x) <- blastNames
    x <- x[!is.na(x$bitscore), ]

    x <- x %>%
      filter(!grepl('Synthetic|clone', stitle, ignore.case = TRUE)) %>%
      arrange(qseqid, desc(bitscore)) %>%
      group_by(qseqid) %>%
      summarize(BID1 = first(staxids),
                BID2 = nth(unique(staxids), 2),
                BID3 = nth(unique(staxids), 3)) %>%
      ungroup()

    # Convert BID1, BID2, BID3 to character type after they are defined
    x <- x %>%
      mutate(
        BID1 = as.character(BID1),
        BID2 = as.character(BID2),
        BID3 = as.character(BID3)
      )

    x <- x %>%
      mutate(BID1 = if_else(str_detect(BID1, ";"), str_split(BID1, ";")[[1]][1], BID1),
             BID2 = if_else(str_detect(BID1, ";"), str_split(BID1, ";")[[1]][2],
                            if_else(str_detect(BID2, ";"), str_split(BID2, ";")[[1]][1], BID2)),
             BID3 = if_else(str_detect(BID2, ";") & length(str_split(BID2, ";")[[1]]) > 1, str_split(BID2, ";")[[1]][2],
                            if_else(str_detect(BID3, ";"), str_split(BID3, ";")[[1]][1], BID3))) %>%
      mutate(BID2 = if_else(is.na(BID2) | BID2 == BID1, NA_character_, BID2),
             BID3 = if_else(is.na(BID3) | BID3 == BID1 | BID3 == BID2, NA_character_, BID3))
    return(x)
  } else {
    return(NULL)
  }
}

allDir <- list.files('./process')

for (sampleDir in allDir) {
  sampleDirPath <- paste0('./process/', sampleDir, '/')
  setwd(sampleDirPath)
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)

  targDir <- './splitSeq10K/'
  filesList <- list.files(targDir)
  targSamples <- paste0(targDir, filesList, '/NtV4_blast.results')
  sampleName <- read.table('sample.name', F)

  blastComb <- foreach(i = targSamples, .combine = rbind, .packages = c('readr', 'dplyr', 'stringr')) %dopar% {
    filterReads(i)
  }

  if (!is.null(blastComb)) {
    blastComb <- blastComb %>%
      rename(Read = qseqid) %>%
      select(Read, BID1, BID2, BID3) %>%
      mutate(Sample = sampleName$V1)

    outFile <- 'blastResTop_3.tsv'
    write.table(blastComb, file = outFile, col.names = T, row.names = F, quote = T, sep = '\t')
  }

  stopCluster(cl)
  print(paste0(sampleDirPath, '________done!'))
  setwd('../../')
}

# Setting the working directory and reading in the data
dir.create('./output/')
setwd("./process/R1_001.fastq.gz")

blast_df <- read.delim("blastResTop_3.tsv", stringsAsFactors = FALSE, sep = "\t")
kraken_df <- read.delim("krakenSelVirReads.tsv", stringsAsFactors = FALSE, sep = ",")
setwd("../../output")

df <- inner_join(blast_df, kraken_df, by = "Read")

# Write the results to a new CSV
write.csv(df, "krakBlastConfReads.csv", row.names = FALSE)
