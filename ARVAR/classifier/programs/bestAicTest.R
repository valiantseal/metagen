library(MASS)
library(pscl)
df = read.csv("288_metaAmpIvar_overlapSnv.csv")
df = read.csv("288_metaAmpIvarCons_overlapSnv.csv")
df = read.csv("Ludy_metaAmpIvar_overlapSnv.csv")
df = read.csv("Ludy_metaAmpIvarNC_overlapSnv.csv")
df = read.csv("Ludy_ampMetaIvarDedup_overlapSnv.csv")
df = read.csv("Ludy_metaAmpIvar_overlapSnv_train.csv")
df = read.csv('Ludy_metaAmpIvar_overlapSnv_RelPos.csv')
df = read.csv("Ludy_metaAmpIvar_overlapSnv_RelPos_RemCont.csv")

colOpt3 = c('ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP')

X <- df[, colOpt3]  # Your independent variable data
y <- factor(df$ConsTest) 



model <- glm(y ~ ., data = X, family = binomial)

# Perform stepwise AIC selection
step_model <- stepAIC(model, direction = "both", scope = list(lower = ~1, upper = ~.), 
                      trace = FALSE)

# Print the selected model
summary(step_model)

###


# Create a data frame of the independent variables

# Create a logistic regression model with all three independent variables
model1 <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV+TOTAL_DP, data = df, family = "binomial")

# Perform a stepwise AIC test to select the best model
step_model  = stepAIC(model1, direction = "backward")
summary(step_model)

step_model  = stepAIC(model1, direction = "both", scope = list(lower = ~1, upper = ~.) , trace = F)
summary(step_model)

# current best AIC ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_QUAL + ALT_RV + TOTAL_DP, family = "binomial", data = df)
# current best AIC ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_QUAL + ALT_RV + TOTAL_DP, family = "binomial", data = df)

step_model <- stepAIC(model1, direction = "both", trace = T)

# Print the selected model
summary(step_model)

# beas AUC model ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV'
## if use on all data without splitting in training and testing 'ALT_FREQ', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'ALT_RV'

model1 <- glm(ConsTest ~ ALT_FREQ+ALT_QUAL+ALT_DP+REF_DP+REF_QUAL+REF_RV+ALT_RV, data = df, family = "binomial")
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
step_model  = stepAIC(model1, direction = "backward")
summary(step_model)

# test current best aic model
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV, data = df, family = "binomial")
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
step_model  = stepAIC(model1, direction = "backward")
summary(step_model)

# previous aic model would be wrong since several variables will be clearly linearly related
# better one would be ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV
# for 288 dataset best model is ALT_FREQ + ALT_QUAL + ALT_DP + REF_RV + ALT_RV

model1 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos+Ref_Al_RelPos, data = metaResDf, family = "binomial")
model1 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos, data = metaResDf, family = "binomial")
summary(model1)

# evaluate with the relative positions
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos + Ref_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)

# try without qualities
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_DP + REF_DP + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model)

# best AIC model: ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)

##
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
model2 =  glm(ConsTest ~ ALT_FREQ + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
model3 = glm(ConsTest ~ ALT_FREQ + ALT_DP + REF_DP +  ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
step_model  = stepAIC(model3, direction = "both" , trace = T, steps = 10000)
summary(step_model)

anova(model2, model1, test = "Chisq")
anova(model3, model1, test = "Chisq")

lrt = pR2(model1)
pR2(model2)
pR2(model3)

sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variables = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}

modRes2 = sumModel(model2)

# try coeficients with best auc including relative position
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + REF_QUAL + REF_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)


# test model with removed contaminated samples
library(MASS)
library(pscl)
df = read.csv("Ludy_metaAmpIvar_overlapSnv_RelPos_RemCont.csv")
dfFilt = df[!(df$Var_Al_RelPos == "NaN"),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))
model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos + Ref_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)
step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ALT_FREQ + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos

best_aic = glm(ConsTest ~ ALT_FREQ + ALT_DP + REF_DP + REF_QUAL + ALT_RV + Var_Al_RelPos, data = dfFilt, family = "binomial")
best_aic_sum = sumModel(best_aic)

write.csv(best_aic_sum , "LogIvarMeta_RelPos_CoefPval_BestAic_RemoveContam.csv", row.names = F)

system("aws s3 cp LogIvarMeta_RelPos_CoefPval_BestAic_RemoveContam.csv s3://abombin/Vivacity/classifier/")
