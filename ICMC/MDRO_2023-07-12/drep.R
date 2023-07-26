rtr_df = read.csv("metadata/rtr_meta_edit.csv")

freq_dat = data.frame(table(rtr_df$Patient_code))

freqSel = freq_dat[(freq_dat$Freq > 1), ]

rtrSel = rtr_df[(rtr_df$Patient_code %in% freqSel$Var1),]

patients = unique(rtrSel$Patient_code)

dir.create("deduplicate_fasta", showWarnings = F)

# need patients id to work correctly

for (patient in patients) {
  curPatient = rtrSel[rtrSel$Patient_code == patient ,]
  bactDirs = unique(curPatient$Bact_id)
  for (bactDir in bactDirs) {
    curBactDir = curPatient[curPatient$Bact_id == bactDir ,]
    if ( nrow(curBactDir) > 1) {
      filesList = list.files(paste0("bactopia/", bactDir, "/output"))
      filesList = filesList[!filesList%in%"bactopia-runs"]
      uuid = unique(curBactDir$uuid)
      selSamples = filesList[grep(paste(uuid, collapse = "|"), filesList)]
      newBactDir = paste0("drep/", bactDir, "/")
      dir.create(newBactDir, recursive = T, showWarnings = F)
      for (curSample in selSamples) {
        inFile = paste0("bactopia/", bactDir, "/output/", curSample, "/main/assembler/", curSample, ".fna.gz")
        outFile = paste0(newBactDir, curSample, ".fna.gz")
        try({
          file.rename(inFile, outFile)
        })

      }
    }
  }
}

runDrep = function(inPath) {
  bactDirs = list.files(inPath, full.names = T)
  for (bactDir in bactDirs) {
    curPatients = list.files(bactDir, full.names =  T)
    for (patient in curPatients) {
      fastaList = list.files(patient)
      if ( length(fastaList) > 1 ) {
        drepDir = paste0(patient, "/drep/")
        system(paste0("rm -rf ", drepDir))
        dir.create(drepDir, recursive = T, showWarnings = F)
        gzip = paste0("gzip -d ", patient, "/*.gz")
        try({
          system(gzip)
        })
        cmd_str = paste0("dRep dereplicate -p 20 -sa 0.995 ", drepDir, " -g ", patient, "/*.fna --checkM_method taxonomy_wf -p 16")
        try({
          system(cmd_str)
        })
      }
    }
  }
  
}

# to run need additionally
# export PATH="/home/ubuntu/extraVol/biosoft/FastANI-master:$PATH"
# checkm data setRoot /home/ubuntu/extraVol/references
runDrep(inPath = "drep")

repId = 
write.table()