# when merging remove time from date and replace unkown and unknown with NA

setwd('/home/ubuntu/ICMC')

## Ahmed
# euhm_cre
samples<-read.table('./Ahmed/euhm_cre/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Babiker'
samples$Project<-'euhm_cre'
samples$Date_processed<-file.info('./Ahmed/euhm_cre/newdir.txt')$ctime
colnames(control2)[1]<-'Sample'
samples<-plyr::join(samples, control2, by="Sample", type='left')
samples$Experiment_Taxa[samples$Experiment_Taxa==""]<-"unknown"
colnames(samples)[5:6]<-c('Gtdb_species', 'Original_species')
samples$Path<-'s3://transfer-files-emory/ICMC/Ahmed/euhm_cre/'
write.csv(samples, 'Babiker_euhm_cre_summary.csv', row.names = F)

#steno-babiker
samples<-read.table('./Ahmed/steno-babiker/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Babiker'
samples$Project<-'steno-babiker'
samples$Date_processed<-file.info('./Ahmed/steno-babiker/newdir.txt')$ctime
colnames(control2)[1]<-'Sample'
samples<-plyr::join(samples, control2, by="Sample", type='left')
samples$Experiment_Taxa[samples$Experiment_Taxa==""]<-"unknown"
colnames(samples)[5:6]<-c('Gtdb_species', 'Original_species')
samples$Path<-'s3://bactopia-outputs/steno-babiker/'
write.csv(samples, 'Babiker_steno-babiker_summary.csv', row.names = F)

# stenotrophomonas-babiker
samples<-read.table('./Ahmed/stenotrophomonas-babiker/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Babiker'
samples$Project<-'stenotrophomonas-babiker'
samples$Date_processed<-file.info('./Ahmed/stenotrophomonas-babiker/newdir.txt')$ctime
setwd('/home/ubuntu/ICMC/Ahmed/stenotrophomonas-babiker')
colnames(mergedId)[1]<-'Sample'
samples<-plyr::join(samples, mergedId, by="Sample", type='left')
colnames(samples)[5]<-c('Gtdb_species')
samples$Original_species<-'Not_provided'
samples$Path<-'Did_not_complete'
write.csv(samples, 'Babiker_stenotrophomonas-babiker_summary.csv', row.names = F)

## Cecile
# cre_witt-cecile
df<-read.csv('/home/ubuntu/ICMC/Cecile/cre_witt-cecile/summary/summary_combined.csv')
samples<-data.frame(unique(df$sample))
colnames(samples)[1]<-'Sample'
samples$Name<-'Cecile'
samples$Project<-'cre_witt-cecile'
samples$Date_processed<-file.info('./Cecile/cre_witt-cecile/bacteria.dirs')$ctime
samples$Gtdb_species<-'Not_tested'
samples$Original_species<-df$organism
samples$Path<-'s3://transfer-files-emory/ICMC/Cecile/cre_witt_results/'
write.csv(samples, 'Cecile_cre_witt-cecile_summary.csv', row.names = F)

# ecoli-cecile
samples<-read.table('./Cecile/ecoli-cecile/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Cecile'
samples$Project<-'ecoli-cecile'
samples$Date_processed<-file.info('./Cecile/ecoli-cecile/newdir.txt')$ctime
samples$Gtdb_species<-'Not_tested'
samples$Original_species<-"Escherichia coli"
samples$Path<-'s3://bactopia-outputs/ICMC/ecoli-cecile/'
write.csv(samples, 'Cecile_ecoli-cecile_summary.csv', row.names = F)

# witt_isolates_09-15-22
samples<-read.table('/home/ubuntu/ICMC/Cecile/witt_isolates_09-15-22/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Cecile'
samples$Project<-'witt_isolates_09-15-22'
samples$Date_processed<-file.info('/home/ubuntu/ICMC/Cecile/witt_isolates_09-15-22/newdir.txt')$ctime
setwd('/home/ubuntu/ICMC/Cecile/witt_isolates_09-15-22')
samples<-cbind(samples, checkedSamples$Experiment_taxa, checkedSamples$Control_taxa)
colnames(samples)[5:6]<-c('Gtdb_species', 'Original_species')
setwd('/home/ubuntu/ICMC')
samples$Path<-'s3://transfer-files-emory/ICMC/Cecile/witt_isolates_09-15-22/'
write.csv(samples, 'Cecile_witt_isolates_09-15-22.csv', row.names = F)

# Alex
# klebastiella_09-02-22
samples<-read.table('/home/ubuntu/ICMC/Page/klebastiella_09-02-22/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Alex'
samples$Project<-'klebastiella_09-02-22'
samples$Date_processed<-file.info('/home/ubuntu/ICMC/Page/klebastiella_09-02-22/newdir.txt')$ctime
setwd('/home/ubuntu/ICMC/Page/klebastiella_09-02-22')
colnames(control2)[1]<-'Sample'
samples<-plyr::join(samples, control2, by="Sample", type='left')
samples$Experiment_Taxa[samples$Experiment_Taxa==""]<-"unknown"
colnames(samples)[5:6]<-c('Gtdb_species', 'Original_species')
setwd('/home/ubuntu/ICMC')
samples$Path<-'s3://transfer-files-emory/ICMC/Page/klebastiella_09-02-22/'
write.csv(samples, 'Alex_klebastiella_09-02-22.csv', row.names = F)

# Sarah
samples<-read.table('/home/ubuntu/ICMC/Sarah/Hflu_sequences/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Sarah'
samples$Project<-'Hflu_sequences'
samples$Date_processed<-file.info('/home/ubuntu/ICMC/Page/klebastiella_09-02-22/newdir.txt')$ctime
samples$Gtdb_species<-'Not_tested'
samples$Original_species<-"Haemophilus influenzae"
samples$Path<-'EC2'
write.csv(samples, 'Sarah_Hflu_sequences_summary.csv', row.names = F)

# Bri
samples<-read.table('/home/ubuntu/ICMC/Bri/newdir.txt')
colnames(samples)[1]<-'Sample'
samples$Name<-'Bri'
samples$Project<-'Bri'
samples$Date_processed<-file.info('/home/ubuntu/ICMC/Bri/newdir.txt')$ctime
setwd('/home/ubuntu/ICMC/Bri/')
colnames(control2)[1]<-'Sample'
samples<-plyr::join(samples, control2, by="Sample", type='left')
samples$Experiment_Taxa[samples$Experiment_Taxa==""]<-"unknown"
colnames(samples)[5:6]<-c('Gtdb_species', 'Original_species')
setwd('/home/ubuntu/ICMC')
samples$Path<-'s3://bactopia-outputs/Bri'
write.csv(samples, 'Bri_summary.csv', row.names = F)
