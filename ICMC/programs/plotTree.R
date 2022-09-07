library(ggtree)
library(ggplot2)

tree<-ape::read.tree('/home/ubuntu/ICMC/Page/klebastiella_09-02-22/bactopia_gtdbtk/pseudomonas-aeruginosa/panOut/bactopia-tools/pangenome/pangenome/myIqtree/boot/core_alignment.fasta.contree')


ggtree::ggtree(tree)+
  geom_tiplab(size=4, align=TRUE, linetype='dashed', linesize=.3)+
  geom_treescale()

#rootTree<-ape::root(tree, outgroup='028_004243', resolve.root = TRUE)


ggtree::ggtree(rootTree)+
  geom_tiplab(size=4, align=TRUE, linetype='dashed', linesize=.3)+
  geom_treescale()

  

ggsave('simpleTree.jpeg', width = 23, height = 16, units = 'in', dpi=300)



ggtree::ggtree(tree, layout="circular")+
  geom_tiplab()+
  geom_treescale()+
  theme(text = element_text(size = 12))



ggsave('simpleTreeCirc.jpeg', width = 14, height = 14, units = 'in', dpi=300)