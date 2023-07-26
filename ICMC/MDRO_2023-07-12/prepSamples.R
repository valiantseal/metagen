library(doParallel)
library(foreach)

df_rtr = read.csv('metadata/metadata_RTR_forbatopia.csv')

df = read.csv("metadata/SRA_IDS_GAMuGSI_updated07052.csv")


renameSpecies<-function(df) {
  
  df$Species = gsub("_", " ", df$Species)
  df$Species = gsub("-", " ", df$Species)
  df$Species_full = NA

  for (i in 1:nrow(df)) {
    if (df$Species[i]=='K. pneumoniae') {
      df$Species_full[i]<-'Klebsiella pneumoniae'
    } else if (df$Species[i]=='P. aeruginosa') {
      df$Species_full[i]<-'Pseudomonas aeruginosa'
    } else if (df$Species[i]=='K. aerogenes') {
      df$Species_full[i]<-'Klebsiella aerogenes'
    } else if (df$Species[i]=='C. freundii') {
      df$Species_full[i]<-'Citrobacter freundii'
    } else if (df$Species[i]=='R. ornithinolytica') {
      df$Species_full[i]<-'Raoultella ornithinolytica'
    } else if (df$Species[i]=='unknown') {
      df$Species_full[i]<-'unknown'
    } else if (df$Species[i]=='P. rettgeri') {
      df$Species_full[i]<-'Providencia rettgeri'
    } else if (df$Species[i]=='E. coli') {
      df$Species_full[i]<-'Escherichia coli'
    } else if (df$Species[i]=='Klebsiella/Enterobacter aerogenes') {
      df$Species_full[i]<-'Klebsiella aerogenes'
    } else {
      df$Species_full[i] = df$Species[i]
    }
  }
  return(df)
}

moveFiles = function(df, inPath) {
  bacteriaList = unique(df$Bact_id)
  for (bacteria in bacteriaList) {
    dfSub = df[(df$Bact_id == bacteria),]
    bacteriaName = unique(dfSub$Species_full)
    bactDir = paste0("bactopia/", bacteria, "/input/")
    dir.create(bactDir, recursive = T, showWarnings = F)
    samplesList = unique(dfSub$uuid)
    for ( curId in samplesList ) {
      curSamplesList = list.files(path = inPath, pattern = curId)
      for (curSample in curSamplesList) {
        combInpath = paste0(inPath, "/", curSample)
        combOutPath = paste0(bactDir, curSample)
        file.rename(combInpath, combOutPath)
        bactTable = paste0("bactopia/", bacteria, "/bacteria.id")
        write.table(bacteriaName, bactTable, col.names = F, row.names = F, quote = T)
      }
    }
  }
}

getFoundSamples = function(df, inPath) {
  df$Found = NA
  for ( i in 1:nrow(df) ) {
    curSample = gsub(" ", "",  df$uuid[i])
    curFilesList = list.files(inPath, pattern = curSample, ignore.case = T)
    curFilesStr = paste(curFilesList, collapse = "__")
    df$Found[i] = curFilesStr
  }
  return(df)
}

compressFastq = function(x, inPath) {
  curFile = paste0(inPath, "/", x)
  outFile = paste0(inPath, "/", x, ".gz")
  cmd_str = paste0('gzip -c ', curFile, " > ", outFile)
  system(cmd_str)
  mv_str = paste0("mv ", curFile, " ", inPath, "/decompressed/")
  system(mv_str)
}

runCompressFastqs = function(inPath) {
  dir.create("original_files/decompressed/", recursive = T, showWarnings = F)
  useCores = 7
  filesList = list.files("original_files", pattern = "\\.fastq$")
  
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)
  
  results<-foreach(i = filesList) %dopar% {
    compressFastq(x = i, inPath = inPath)
  }
  
  parallel::stopCluster(cl = cl)
}

# transfer samples from metadata, needs testing
moveMetaFiles = function(df, inPath) {
  for ( i in 1:nrow(df) ) {
    curBactDir = df$Bact_id[i]
    curBactName = df$Species_full[i]
    bactDirIn = paste0("bactopia/", curBactDir, "/input/")
    curFiles = strsplit(df$Found[i], "__")[[1]]
    for (curFile in curFiles) {
      inName = paste0("original_files/", curFile)
      if (grepl(".fastq.gz", curFile)) {
        if (!file.exists(paste0("bactopia/", curBactDir))) {
          dir.create(bactDirIn, recursive = T, showWarnings = F)
          bactTable = paste0("bactopia/", curBactDir, "/bacteria.id")
          write.table(curBactName, bactTable, col.names = F, row.names = F, quote = T)
        }
        combInpath = paste0(inPath, "/", curFile)
        combOutPath = paste0(bactDirIn, curFile)
        file.rename(combInpath, combOutPath)
      } else if (grepl(".fasta", curFile)) {
        fastaInput = paste0("bactopia_fasta/", curBactDir, "/input/")
        fastaBactDir = paste0("bactopia_fasta/", curBactDir)
        dir.create(fastaInput, recursive = T, showWarnings = F)
        combInpath = paste0(inPath, "/", curFile)
        combOutPath = paste0(fastaInput, curFile)
        file.rename(combInpath, combOutPath)
        bactTable = paste0(fastaBactDir, "/bacteria.id")
        write.table(curBactName, bactTable, col.names = F, row.names = F, quote = T)
      }
    }
  }
}


# process sra files
df = renameSpecies(df)
df$Bact_id = tolower(gsub(" ", "_", df$Species_full))
df$uuid = df$SRA.Number
moveFiles(df = df, inPath = "input")

# process files custom
df_rtr = renameSpecies(df_rtr)
df_rtr$Bact_id = tolower(gsub(" ", "_", df_rtr$Species_full))
df_rtr$uuid = df_rtr$genome_ID
#write.csv(df_rtr, "metadata/rtr_meta_edit.csv", row.names = F)
df_rtr = getFoundSamples(df = df_rtr, inPath = "original_files/")

df_rtr_missing = df_rtr[(df_rtr$Found == "") | (is.na(df_rtr$Found)),]
write.csv(df_rtr_missing, 'metadata_RTR_forbatopia_missing.csv', row.names = F)
# system("aws s3 cp metadata_RTR_forbatopia_missing.csv s3://transfer-files-emory/ICMC/Ahmed/MDRO_2023-07-12/")

# compress fastqs
runCompressFastqs(inPath = "original_files")
# move rtr seqs
moveMetaFiles(df = df_rtr, inPath = "original_files")

# check if new samples have fastqs for missing and for fasta files
missDf = read.csv("metadata_RTR_forbatopia_missing.csv")
missDf = getFoundSamples(df = missDf, inPath = "original_files_2/")

fasta = getFoundSamples(df = df_rtr, inPath = "original_files_2/")