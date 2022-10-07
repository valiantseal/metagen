ca1Stress<-allMarkRenClust0.3ResGroups_2022.10.05[(allMarkRenClust0.3ResGroups_2022.10.05$Cell_Type=='CA3' & allMarkRenClust0.3ResGroups_2022.10.05$cluster=='Stress'),]

genesStress<-unique(ca1Stress$gene)

ca1Control<-allMarkRenClust0.3ResGroups_2022.10.05[(allMarkRenClust0.3ResGroups_2022.10.05$Cell_Type=='CA3' & allMarkRenClust0.3ResGroups_2022.10.05$cluster=='Control'),]

genesControl<-unique(ca1Control$gene)

identical(genesControl, genesStress)

pStress<-unique(ca1Stress$p_val)
pControl<-unique(ca1Control$p_val)

identical(pStress, pControl)

w1<-data.frame(genesStress, genesControl)
w2<-w1[!(w1$genesStress==w1$genesControl),]


q1<-data.frame(cbind(pStress, pControl))

q2<-q1[!(q1$pStress == q1$pControl),]

ca1Control$fdr<-p.adjust(ca1Control$p_val, 'bonferroni')


for (i in clusters ){
  ca1Stress<-allMarkRenClust0.3ResGroups_2022.10.05[(allMarkRenClust0.3ResGroups_2022.10.05$Cell_Type==i & allMarkRenClust0.3ResGroups_2022.10.05$cluster=='Stress'),]
  ca1Control<-allMarkRenClust0.3ResGroups_2022.10.05[(allMarkRenClust0.3ResGroups_2022.10.05$Cell_Type==i & allMarkRenClust0.3ResGroups_2022.10.05$cluster=='Control'),]
  genesStress<-unique(ca1Stress$gene)
  genesControl<-unique(ca1Control$gene)
  print(identical(genesControl, genesStress))
}

r1<-RNA.combined.norm@assays$RNA@counts
r2<-data.frame(rowSums(r1))