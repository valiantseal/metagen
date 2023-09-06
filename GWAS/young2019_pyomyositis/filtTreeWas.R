library(treeWAS)

twOut = readRDS("/home/ubuntu/extraVol/GWAS/young2019_pyomyositis/treeWAS/tr_output.rds")

# terminal = data.frame(twOut[["terminal"]][["p.vals"]])
terminal = data.frame(twOut$terminal$p.vals)
colnames(terminal)[1] = "p.val"
terminal$AL = row.names(terminal)

simmult = data.frame(twOut$simultaneous$p.vals)
colnames(simmult)[1] = "p.val"
simmult$AL = row.names(simmult)

subseq = data.frame(twOut$subsequent$p.vals)
colnames(subseq)[1] = "p.val"
subseq$AL = row.names(subseq)

adjPFilter = function(df) {
  dfFilt = df
  dfFilt$Bonfer_p = p.adjust(dfFilt$p.val, method = 'bonferroni')
  dfSign = dfFilt[dfFilt$Bonfer_p < 0.05,]
  rownames(dfSign) = NULL
  return(dfSign)
}

termSign = adjPFilter(terminal)
simmultSign = adjPFilter(simmult)
subseqSign = adjPFilter(subseq)