mySnv = read.csv("combDatLudyTrueSnvs.csv")

vivSnv = read.csv("/home/ubuntu/extraVol/ARVAR/Vivacity/Ludy_metaseq/combDatFilt.csv")

mySamples = unique(mySnv$Sample)
vivSamples = unique(vivSnv$Sample)

length(vivSamples[vivSamples%in%mySamples])

mySnvFilt = mySnv[(mySnv$Sample%in%vivSamples),]

length(unique(mySnv$Samp_Pos_Ref_Alt))
length(unique(mySnvFilt$Samp_Pos_Ref_Alt))
length(unique(vivSnv$Samp_Pos_Ref_Alt))


mySnvs = unique(mySnv$Samp_Pos_Ref_Alt)
mySnvsFilt = unique(mySnvFilt$Samp_Pos_Ref_Alt)
vivSnvs = unique(vivSnv$Samp_Pos_Ref_Alt)

length(mySnvs[mySnvs%in%vivSnvs])
length(mySnvsFilt[mySnvsFilt%in%vivSnvs])



length(mySnvs[mySnvs%in%vivSnvs]) / length(unique(mySnv$Samp_Pos_Ref_Alt)) * 100

length(mySnvs[mySnvs%in%vivSnvs]) / length(unique(mySnvFilt$Samp_Pos_Ref_Alt)) * 100

length(mySnvs[mySnvs%in%vivSnvs]) / length(unique(vivSnv$Samp_Pos_Ref_Alt)) * 100