library(ggvenn)
library(ggpubr)

mySnv = read.csv("combDatLudyTrueSnvs.csv")

vivSnv = read.csv("/home/ubuntu/extraVol/ARVAR/Vivacity/Ludy_metaseq/combDatFilt.csv")

mySamples = unique(mySnv$Sample)
vivSamples = unique(vivSnv$Sample)
length(mySamples)
length(vivSamples)

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


unique(mySnv$Sample[!(mySnv$Sample%in%vivSnv$Sample)])


# make venn Diagrams
genes = list("Vivacity_SNVs" = unique(vivSnv$Samp_Pos_Ref_Alt), "Ivar/LogClass_SNVs" = unique(mySnv$Samp_Pos_Ref_Alt))
vennPlot <-ggvenn(genes, text_size=5, set_name_size=8)+
  theme(text = element_text(size = 26)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -6))


png(file = "Venn_LudyAllMetaseq_Vivacity_LogClass.png", width = 9, height = 9, units = 'in', res = 300)
print(vennPlot)
dev.off()