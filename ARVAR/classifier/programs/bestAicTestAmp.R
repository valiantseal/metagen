library(MASS)
library(pscl)

df = read.csv("result_tables/Ludy_ampDedupMetaIvar_overlapSnv_RelPos.csv")

dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
colnames(dfFilt)

completeModel <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV +Var_Al_RelPos+Ref_Al_RelPos , data = dfFilt, family = "binomial")
summary(completeModel)
completeStep  = stepAIC(completeModel, direction = "both" , trace = T, steps = 10000)
summary(completeStep) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos

model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos + Ref_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos

# test standard parameter without position
model1 <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV, data = dfFilt, family = "binomial")
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV

# make test based on the previous results and add position
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV + Var_Al_RelPos + Ref_Al_RelPos, data = dfFilt, family = "binomial")
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos

# testing variables separately or all together gives the same results
anova(model1, completeModel, test = "Chisq")

# coeficients with the best auc score
aucModel = glm(ConsTest ~ ALT_FREQ+ALT_DP+REF_QUAL+ALT_RV+Ref_Al_RelPos , data = dfFilt, family = "binomial")
summary(aucModel)

# get coeficients table
orig_aic <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV, data = dfFilt, family = "binomial")
pos_aic = glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos, data = dfFilt, family = "binomial")

sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variables = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}

orig_aic_sum = sumModel(orig_aic)

pos_aic_sum = sumModel(pos_aic)

write.csv(orig_aic_sum, "ampSeq_Dedup_Coef_BestOriginalAic.csv", row.names = F)
write.csv(pos_aic_sum, "ampSeq_Dedup_Coef_BestRelPosAic.csv", row.names = F)

system("aws s3 cp ampSeq_Dedup_Coef_BestOriginalAic.csv s3://abombin/Vivacity/classifier/")
system("aws s3 cp ampSeq_Dedup_Coef_BestRelPosAic.csv s3://abombin/Vivacity/classifier/")


model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)



# removed contamination 
df = read.csv("result_tables/Ludy_ampDedupMetaIvar_overlapSnv_RelPosRemCont.csv")

dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
colnames(dfFilt)

completeModel <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV +Var_Al_RelPos+Ref_Al_RelPos , data = dfFilt, family = "binomial")
summary(completeModel)
completeStep  = stepAIC(completeModel, direction = "both" , trace = T, steps = 10000)
summary(completeStep) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + Ref_Al_RelPos

# original without position
model1 <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV, data = dfFilt, family = "binomial")
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV