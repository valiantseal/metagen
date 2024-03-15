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
  if (!file.exists(df_path)) {
    warning(paste("File does not exist:", df_path))
    return(NULL) 
  }

  x <- read_delim(df_path, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)

  if (nrow(x) == 0) {
    warning(paste("Dataframe is empty after loading file:", df_path))
    return(NULL)
  }

  colnames(x) <- blastNames
  x <- x[!is.na(x$bitscore), ]

  if (nrow(x) == 0) {
    warning(paste("Dataframe is empty after subsetting by bitscore:", df_path))
    return(NULL)
  }

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

  return(x)
}

registerDoParallel(cores = useCores) 

allDir <- list.files('./process', full.names = TRUE)

for (sampleDirPath in allDir) {
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)

  targDir <- paste0(sampleDirPath, '/splitSeq10K/')
  filesList <- list.files(targDir, pattern = "*.reads_dir", full.names = TRUE)
  targSamples <- sapply(filesList, function(x) paste0(x, '/NtV4_blast.results'), simplify = FALSE)
  targSamples <- unlist(targSamples)
  sampleName <- read.table(paste0(sampleDirPath, '/sample.name'), col.names = c("V1"))

  existingTargSamples <- Filter(file.exists, targSamples)

  results <- foreach(i = existingTargSamples, .combine = rbind, .packages = c('readr', 'dplyr', 'stringr', 'tidyr')) %dopar% {
    tryCatch({
      filterReads(i, blastNames)
    }, error = function(e) {
      warning(paste("Error processing file:", i, "\nError message:", e$message))
      NULL 
    })
  }

  if (!is.null(results)) {
    results <- results %>%
      rename(Read = qseqid) %>%
      select(Read, BID1, BID2, BID3) %>%
      mutate(Sample = sampleName$V1)

    outFile <- paste0(sampleDirPath, '/blastResTop_3.tsv')
    write.table(results, file = outFile, col.names = TRUE, row.names = FALSE, quote = TRUE, sep = '\t')
  }

  stopCluster(cl)
  print(paste0(sampleDirPath, '________done!'))
}
