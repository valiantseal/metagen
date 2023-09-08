pheno = read.delim("pyseer_input/bin_phenotype.tsv")
colnames(pheno)[1] = "Name"
write.csv(pheno, "bin_phenotype.csv", row.names = F)

dir.create("Scoary", showWarnings = F)

geneMat = "pangenome_panaroo/bactopia-runs/pangenome-20230826-101558/panaroo/gene_presence_absence_roary.csv"

cmd_str = "scoary -g pangenome_panaroo/bactopia-runs/pangenome-20230826-101558/panaroo/gene_presence_absence_roary.csv \\
-t bin_phenotype.csv \\
 --threads 30 \\
 -o Scoary"

#system(cmd_str)

cmd_str1 = "scoary -g pangenome_panaroo/bactopia-runs/pangenome-20230826-101558/panaroo/gene_presence_absence_roary.csv \\
-t bin_phenotype.csv \\
 --threads 30 \\
 -s 2 \\
 -o Scoary"

system(cmd_str1)

# pvalues for cmd and cmd1 are identical