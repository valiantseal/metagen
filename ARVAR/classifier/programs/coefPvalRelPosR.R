df = read.csv('Ludy_metaAmpIvar_overlapSnv_RelPos.csv')

dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
model2 =  glm(ConsTest ~ ALT_FREQ + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
model3 = glm(ConsTest ~ ALT_FREQ + ALT_DP + REF_DP +  ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")

sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variables = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}

modRes2 = sumModel(model2)
modRes3 = sumModel(model3)

write.csv(modRes2, "LogIvarMeta_RelPos_CoefPval_SignVars.csv", row.names = F)

write.csv(modRes3, "LogIvarMeta_RelPos_CoefPval_ExclQual.csv", row.names = F)

system("aws s3 cp LogIvarMeta_RelPos_CoefPval_SignVars.csv s3://abombin/Vivacity/classifier/")
system("aws s3 cp LogIvarMeta_RelPos_CoefPval_ExclQual.csv s3://abombin/Vivacity/classifier/")
