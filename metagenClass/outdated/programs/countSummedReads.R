library(dplyr)

df0<- read.table('/home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/blastNtSummary/B2E22-002A_S1_L001.txt', F)

dfSum0<- df0 %>% group_by(V2) %>%
  summarise(Freq = sum(V1))

df1<- read.table('/home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/blastNtSummary/B2E22-003A_S2_L001.txt', F)

dfSum1<- df1 %>% group_by(V2) %>%
  summarise(Freq = sum(V1))


df2<- read.table('/home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/blastNtSummary/B2E22-004A_S3_L001.txt', F)

dfSum2<- df2 %>% group_by(V2) %>%
  summarise(Freq = sum(V1))

