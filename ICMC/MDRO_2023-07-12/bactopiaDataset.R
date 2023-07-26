checkBactopiaDb = function(inPath) {
  curData = list.files(inPath)
  curData = gsub("_", "-", curData)
  bactopData = list.files('/home/ubuntu/bactopia/datasets/species-specific')
  
  missData = curData[!curData%in%bactopData]
  if (length(missData > 0)) {
    for (bactDir in missData) {
      curBactDir = gsub("-", "_", bactDir)
      curBacteria = read.table(paste0(inPath, "/", curBactDir, "/bacteria.id"))
      curBacteria = as.character(curBacteria$V1)
      cmd_str = paste0("bactopia datasets --species '", as.character(curBacteria), "' --outdir /home/ubuntu/bactopia/datasets --cpus 8")
      try({
        system(cmd_str)
      })
    }
  }
}

finalCheck = function(inPath) {
  curData = list.files(inPath)
  bactopData = list.files('/home/ubuntu/bactopia/datasets/species-specific')
  curData = gsub("_", "-", curData)
  missData = curData[!curData%in%bactopData]
  if (length(missData > 0)) {
    for (bactDir in missData) {
      print(paste0(bactDir, " is missing"))
    } 
  } else {
    print("All datasets are present")
  }
}

checkBactopiaDb(inPath = "bactopia")
finalCheck(inPath = "bactopia")

checkBactopiaDb(inPath = "bactopia_fasta")
finalCheck(inPath = "bactopia_fasta")