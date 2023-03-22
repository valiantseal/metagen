assocRes = read.delim('./gwas_output/simple_assoc_pyseer.txt')
lmmRes = read.delim('./gwas_output/lmm_unitig_pyseer.txt')

print(nrow(assocRes))
assocRes$bonf_p = p.adjust(assocRes$lrt.pvalue)
print(length(assocRes$bonf_p[assocRes$bonf_p < 0.05]))

print(nrow(lmmRes))
lmmRes$bonf_p = p.adjust(lmmRes$lrt.pvalue)
print(length(lmmRes$bonf_p[lmmRes$bonf_p < 0.05]))

