primers<-read.delim('../../references/DENV1_reference.bed', F, skip = 1)

fasta <- phylotools::read.fasta("../../references/DENV1_reference.fasta")

fasta$seq.name

fasta$seq.name <- gsub("\\ .*", "", fasta$seq.name)

identical(unique(primers$V1), fasta$seq.name)


write.table(primers, '../../references/DENV1_reference_edit.bed', row.names = F, col.names = F, quote = F, sep = '\t')

phylotools::dat2fasta(fasta, "../../references/DENV1_reference_edit.fasta")
