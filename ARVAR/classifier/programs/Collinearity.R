library(MASS)
library(pscl)
library(car)

# current ivar metaseq model
df = read.csv("result_tables/Ludy_metaAmpIvar_overlapSnv_RelPos_RemCont.csv")
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt = dfFilt[!(is.na(dfFilt$Var_Al_RelPos)),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

best_aic = glm(ConsTest ~ ALT_FREQ + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
summary(best_aic)
vif_values <- vif(best_aic)
print(vif_values)

# try with all variables from 
allMetaseq = glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + REF_RV + TOTAL_DP +  Var_Al_RelPos + Ref_Al_RelPos, data = dfFilt, family = "binomial")
print(vif(allMetaseq))

step_model  = stepAIC(model1, direction = "both" , trace = F, steps = 10000)
summary(step_model)

# current ivar ampseq
df = read.csv("result_tables/Ludy_ampDedupMetaIvar_overlapSnv_RelPosRemCont.csv")
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt = dfFilt[!(is.na(dfFilt$Var_Al_RelPos)),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
completeModel <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV+Var_Al_RelPos+Ref_Al_RelPos , data = dfFilt, family = "binomial")
summary(completeModel)
completeStep  = stepAIC(completeModel, direction = "both" , trace = F, steps = 10000)
summary(completeStep) 

best_aic1 = glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos , data = dfFilt, family = "binomial")
best_aic2 = glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV + Var_Al_RelPos , data = dfFilt, family = "binomial")

print(vif(best_aic1))
print(vif(best_aic2))
print(vif(completeModel))

# current lofreq metaseq
df = read.csv("result_tables/Ludy_metaAmp_ViralrecLofreq_overlapSnv_RemCont.csv")
dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos, data = dfFilt, family = "binomial")
print(vif(multivarModel))
summary(multivarModel)
