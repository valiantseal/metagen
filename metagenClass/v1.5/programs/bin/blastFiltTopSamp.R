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
    return(NULL)  # Skip this file if it doesn't exist
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
      mutate(BID1 = if_else(str_detect(BID1, ";"),
                            str_extract(BID1, "^[^;]+"),
                            BID1),
             BID2 = if_else(BID2 == BID1, NA_character_, BID2),
             BID3 = if_else(BID3 == BID1 | BID3 == BID2, NA_character_, BID3))
  
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
