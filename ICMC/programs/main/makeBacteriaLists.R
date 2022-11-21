# read metadata
metadata<-readr::read_csv("./metadata/metadata.csv")

# format metadata
processMetadata<-function(x){
  metadatID<-metadata[, c('SRA','Species_full')]
  colnames(metadatID)<-c("uuid", "Control_taxa")
  metadatID$Control_taxa[metadatID$Control_taxa==""]<-'unknown'
  metadatID$Control_taxa[metadatID$Control_taxa=="Clostridium difficile"]<-"Clostridioides difficile"
  return(metadatID)
}

metadataID<-processMetadata(metadata)

#extract matched

combineMatched<-unique(metadataID)

# write summary for app
idSummary<-metadataID
colnames(idSummary)[2]<-'Original_species'
idSummary$Gtdb_species<-'Not_tested'
write.csv(idSummary, 'idSummary.csv', row.names = F)

colnames(combineMatched)[2]<-'Final_taxa'

taxa<-unique(combineMatched$Final_taxa)

# make directories for output

dir.create("./test_gtdbtk/buckets/", recursive = T)
dir.create("./test_gtdbtk/taxa/")

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