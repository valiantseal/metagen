df<-read.delim('MLS_teammates_metaEdit.tsv', T, sep='\t')


df$date<-as.Date(df$date, format="%m/%d/%y")
df$date<-format(df$date, format="%Y-%m-%d")

df$date_submitted<-df$date

write.table(df, 'MLS_teammates_metaEdit.tsv', row.names = F, quote = F, sep = '\t')


df1<-df[!is.na(df$date),]


write.table(df1, 'AY4all_metaEdit.tsv', row.names = F, quote = F, sep = '\t')