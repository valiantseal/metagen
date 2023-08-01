system("rm -rf ./panOut")
system("bactopia --wf pangenome --bactopia ./bactopia_output --outdir ./panOut --max_cpus 16 --skip_phylogeny --skip_recombination --run_name 'klebsiella_pneumoniae' ")

curDirs = list.files(paste0("panOut/bactopia-runs"), full.names = T)
# change next line when run with normal bactopia runs
for (curDir in curDirs) {
  setwd(curDir)
  iqCommand<-paste0('iqtree -s core-genome.aln.gz -m TEST -nt 16 -bb 1000')
  system(iqCommand)
  setwd("../../../")
}

