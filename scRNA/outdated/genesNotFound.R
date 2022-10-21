q1<-RNA.combined.norm@assays$RNA@data@Dimnames
all_genes<-q1[[1]]
q2<-q1[[2]]

pyr_markers<-read.csv('../data/pyr_markers.csv')

pyr<-pyr_markers$markers

all_select_genes<-c("Rbfox3", "Snap25", "Syt1", "Slc17a6", "Grin2b", "Gad1", "Gad2","Slc1a3", "Gfap", "Aldoc", 
                    "Glul", "Rarres2", "Slc6a13", "Plp1", "Vcan", "Hmha1", "Iba-1","Alf-1","CX3CR1", "P2ry12", "Reln", "Ndnf")

all_neuron_subtypes<-c("Reln", "Ndnf", "Mpped1", "Satb2", "Map3k15", "Gram4","Cdh24", "Npy2r",  "Prox1", "Glis3")

additional_markers<-c("Cx3cr1", "Ephb1", "Vwc2I", "Csf2rb2", "Fibcd1", "Wfs1","DCN",  "Slc4a4", "Plk5" , "Cntn6", "Amigo2")
addMarkSet2<-c("Kcnh5", "Bok", "CSmd3","PCDH9", "Akt2", "Acsl3", "Trf", "Ptgds", "Hexb", "Kcnc1", "Npy1r", "Calb2", "Calb1")
addMarkSet3<-c("Olig2", "Gramd4", "Vwc2l")

genes<-c('Rbfox3', 'Snap25', 'Syt1', 'Fibcd1', 'Mpped1', 'Satb2', 'Amigo2', 'Cntn6', 'Kcnh5', 'Ephb1', 'Glis3', 'Plk5', 'Slc4a4', 'Gad1', 'Gad2', 'Plp1', 'Ptgds', 'Trf', 'Vcan', 
         'Cx3cl1', 'Hexb', 'Ptry12', 'Ndnf', 'Reln')


all_used_markers<-unique(c(pyr, all_select_genes, all_neuron_subtypes, additional_markers, addMarkSet2, addMarkSet3, genes))

MarkersNotFound<-data.frame(all_used_markers[!(all_used_markers%in%all_genes)])
colnames(MarkersNotFound)<-'Not_found_genes'


write.csv(MarkersNotFound, 'Not_found_markers.csv', row.names = F)