library(readr)

prepGtdb = function(bactDir) {
  setwd(bactDir)
  dir.create("gtdbtk_input", showWarnings = F)
  samplesList = list.files("bactopia_output")
  for (curSample in samplesList) {
    inFile = paste0('bactopia_output/', curSample, "/main/assembler/", curSample, ".fna.gz")
    outFile = paste0("gtdbtk_input/", curSample, ".fna.gz")
    file.copy(inFile, outFile)
  }
}

# needs conda activate gtdbtk-2.1.0
runGtdb = function() {
  dir.create('gtdbtk_output', showWarnings = F)
  cmd_str = "gtdbtk classify_wf --genome_dir gtdbtk_input --out_dir gtdbtk_output --extension gz --cpus 20 --skip_ani_screen"
  system(cmd_str)
}

#prepGtdb(bactDir = "pangenome/klebsiella_pneumoniae/")

runGtdb()


gtdb<-readr::read_delim("gtdbtk_output/gtdbtk.bac120.summary.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# separate main gtdbtk classification
processGtdb<-function(x){
  sepNames<-data.frame(do.call('rbind', strsplit(as.character(x$classification),';',fixed=TRUE)))
  gtdbName<-data.frame(cbind(x$user_genome, sepNames$X7, sepNames$X6))
  colnames(gtdbName)<-c("uuid", "Experiment_taxa", "Experiment_genus")
  gtdbName$Experiment_taxa<-gsub("s__", "", gtdbName$Experiment_taxa)
  gtdbName$Experiment_taxa[gtdbName$Experiment_taxa==""]<-'unknown'
  return(gtdbName)
}

# get main species from gtdb classification
gtdb_result<-processGtdb(gtdb)
gtd_exclude = gtdb_result[!gtdb_result$Experiment_taxa == "Klebsiella pneumoniae",]
excludeList = gsub('.fna', "", gtd_exclude$uuid)

moveSamples = function(excludeList) {
  for ( curSample in excludeList) {
    inFile = paste0('bactopia_output/', curSample, "/cd ")
    cmd_str = paste0("mv ", inFile, " exclude_samples/")
    print(cmd_str)
    system(cmd_str)
  }
}


moveSamples(excludeList = excludeList)