setwd("./")

# read gtdbtk data
gtdb<-readr::read_delim("gtdbtk/output/gtdbtk.bac120.summary.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

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

# read metadata
metadata<-readr::read_csv("./metadata/metadata.csv")

# format metadata
processMetadata<-function(x){
  metadatID<-metadata[, c('uuid','Species_full')]
  colnames(metadatID)<-c("uuid", "Control_taxa")
  metadatID$Control_taxa[metadatID$Control_taxa==""]<-'unknown'
  metadatID$Control_taxa[metadatID$Control_taxa=="Clostridium difficile"]<-"Clostridioides difficile"
  return(metadatID)
}

metadataID<-processMetadata(metadata)

# merge metadata and gtdb
mergedData<-plyr::join(gtdb_result, metadataID, by="uuid", type="left", match="first")
write.csv(mergedData, './run_info/taxaInfo.csv', row.names = F)

# check if gtdb or original taxa are not identified
checkMatch<-function(x){
  df<-x
  df$Final_taxa<-'unknown'
  for (i in 1:nrow(df)){
    if (df$Experiment_taxa[i]==df$Control_taxa[i]){
      df$Final_taxa[i]=df$Experiment_taxa[i]
    } else if (df$Experiment_taxa[i]=='unknown') {
      df$Final_taxa[i]=df$Control_taxa[i]
    } else if (df$Control_taxa[i]=='unknown') {
      df$Final_taxa[i]=df$Experiment_taxa[i]
    }
  }
  firstMatch<-df
}

checkedSamples<-checkMatch(mergedData)

# extract matched samples
firstMatch<-checkedSamples[!(checkedSamples$Final_taxa=='unknown'),]
# extract mismatch
firstMismatch<-mergedData[!(mergedData$uuid%in%firstMatch$uuid),]

# check spelling
checkSpelling<-function(x){
  spellCheck<-x
  spellCheck$LS<-RecordLinkage::levenshteinSim(spellCheck$Experiment_taxa, spellCheck$Control_taxa)
  spellCheck$Final_taxa <- ifelse((spellCheck$LS>0.75), spellCheck$Control_taxa, "unknown")
  checkPass<-spellCheck[!(spellCheck$Final_taxa=="unknown"),]
  matched_samples=checkPass[, c('uuid', 'Experiment_taxa', 'Experiment_genus', 'Control_taxa', 'Final_taxa')]
  matched_samples<-matched_samples[!(is.na(matched_samples$uuid)),]
  if (nrow(matched_samples>0)){
    rownames(matched_samples)<-1:nrow(matched_samples)
  }
  return(matched_samples)
}

#extract matched
secondMatch<-checkSpelling(firstMismatch)

matched_1_2<-rbind(firstMatch, secondMatch)

# extract mismatch
secondMismatch<-mergedData[!(mergedData$uuid%in%matched_1_2$uuid),]

# check if Gtdb names corresponds to ncbi names
checkNCBI<-function(x){
  ncbi<-read.delim('/home/ubuntu/gtdbtk/conversion/dictionary/gtdbtkNcbiDictSpecies.tsv', T, sep='\t')
  colnames(ncbi)[1]<-'Experiment_taxa'
  ncbiMatch<-plyr::join(x, ncbi, by='Experiment_taxa', match='first', type='left')
  names(ncbiMatch)[names(ncbiMatch) == "ncbi_organism_name"] <- "Final_taxa"
  thrirdMatch<-ncbiMatch[, c("uuid","Experiment_taxa","Experiment_genus","Control_taxa","Final_taxa" )]
  return(thrirdMatch)
}

#extract matched
thirdMatch<-checkNCBI(secondMismatch)

combineMatched<-rbind(firstMatch, secondMatch, thirdMatch)

taxa<-unique(combineMatched$Final_taxa)

# write files
writeSamples<-function(){
  for (taxon in taxa){
    df<-combineMatched[(combineMatched$Final_taxa==taxon),]
    samples<-data.frame(df$uuid)
    bucket<-tolower(gsub(" ", "-", taxon))
    outdir_samples<-paste0("./test_gtdbtk/buckets/",bucket)
    outdir_taxa<-paste0("./test_gtdbtk/taxa/",bucket)
    write.table(samples, file= outdir_samples, col.names = F, row.names = F, quote = F)
    write.table(taxon, file= outdir_taxa, col.names = F, row.names = F, quote = F)
  }
}

writeSamples()

# plan for a new program
# match metadata and gtdb
# if metadata is empty or unkown replace metadata with gtdb and add to log
# if gtdb is empty or unkown replace with metadata and add to log
# if metadata == gtdb pass
# for mismatch check if it is missplelling
# if misspelling correct and pass

