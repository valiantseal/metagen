df <- read.csv("Scoary/Phenotype_08_09_2023_1046.results.csv")

dfSel = df[df$Bonferroni_p < 0.05,]

dfSel = dfSel[!dfSel$Annotation == "hypothetical protein",]