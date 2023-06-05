mySnv = read.csv("combDatLudyTrueSnvs.csv")

# my Samples
mySumTab = data.frame(table(mySnv$ExactSamp))
median(mySumTab$Freq)
subSamp = mySumTab[(mySumTab$Freq >= 46) & (mySumTab$Freq <= 46),]
targSamp = as.character(head(subSamp$Var1, 10))

#

testDat = read.csv('/home/ubuntu/extraVol/ARVAR/Vivacity/Ludy_metaseq/combDatFilt.csv')
testSumTab = data.frame(table(testDat$ExactSamp))
median(testSumTab$Freq)

subTest= testSumTab[(testSumTab$Freq >= 43) & (testSumTab$Freq <= 43),]


selSamp = subTest[(subTest$Var1%in%subSamp$Var1),]

selSamp1 = subTest[(subTest$Var1%in%targSamp),]

# transfer files 

# ivarClass
for (curSamp in targSamp) {
  filesList = list.files("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/output/variants/bowtie2", full.names = T, pattern = curSamp)
  for (curFile in filesList) {
    cmd_str = paste0("aws s3 cp ", curFile, " s3://abombin/Vivacity/ManCheckSamp/LogClass/bams/")
    system(cmd_str)
  }
}

selDat = mySnv[(mySnv$ExactSamp%in%targSamp),]

write.csv(selDat, "LogClassSNV_forManCheck.csv", row.names = F)
system("aws s3 cp LogClassSNV_forManCheck.csv s3://abombin/Vivacity/ManCheckSamp/LogClass/")

# lofreq output

for (curSamp in targSamp) {
  filesList = list.files("/home/ubuntu/extraVol/ARVAR/Vivacity/Ludy_metaseq/process", full.names = T, pattern = curSamp)
  for (curFile in filesList) {
    cmd_str = paste0("aws s3 cp --recursive ", curFile, " s3://abombin/Vivacity/ManCheckSamp/LofreqFilt/output/", curSamp, "/")
    system(cmd_str)
  }
}

selDat = testDat[(testDat$ExactSamp%in%targSamp),]
write.csv(selDat, "LofreqFiltSNV_forManCheck.csv", row.names = F)
system("aws s3 cp LogClassSNV_forManCheck.csv s3://abombin/Vivacity/ManCheckSamp/LofreqFilt/")