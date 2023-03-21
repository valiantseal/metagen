tree<-ape::read.tree('iqtree/core.contree')

tree$tip.label<-gsub('_', '-', tree$tip.label)


ape::write.tree(tree, 'iqtree/relab_core.contree')