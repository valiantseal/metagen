library(ggtree)
library(ggplot2)
library(ggbreak)


setwd('C:/Users/abomb/OneDrive - Emory University/virus/data/patient288')

metaDat<-read.delim('nextstrain__metadata.tsv', T)

metaDat$location[1]<-'Wuhan'

metaDat$Tips<-''

for ( i in 1:nrow(metaDat)){
  if (metaDat$strain[i] == 'USA/GA-EHC-2884X/2021') {
    metaDat$Tips[i] = 'Day 0'
  } else if (metaDat$strain[i] == 'USA/GA-EHC-2885Y/2021') {
      metaDat$Tips[i] = 'Day 31'
  } else if (metaDat$strain[i] == 'USA/GA-EHC-2886Z/2021') {
    metaDat$Tips[i] = 'Day 44'
    }
}

metaDat$location <- factor(metaDat$location, levels = c("patient288", "Georgia"))

group.colors <- c(patient288 = "#D55E00", Georgia = "#0072B2")




tree<-ape::read.tree('nextstrain__tree.nwk')

plotTree<-ggtree(tree)

metaPlot<-plotTree%<+% metaDat + 
  geom_tippoint(aes(color = location), size=2)+
  scale_x_break(c(1, 20)) +
  scale_x_continuous(breaks=seq(0, 45, 5))+
  theme_tree2()+
  scale_color_manual(values=group.colors)+
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())+
  xlab("Mutations")+
  theme(text = element_text(size = 18))+
  theme(legend.position = "none")+
  geom_tiplab(aes(label = Tips), size=3)

ggsave(filename = 'Patient_288_tree.tiff', plot = metaPlot, width = 12, height = 8, units = 'in', dpi = 600, device = 'tiff')
