library(MASS)
library(pROC)
library(caret)
df = read.csv("result_tables/Ludy_metaAmp_ViralrecLofreq_overlapSnv_RemCont.csv")

dfFilt = df[!is.na(df$Var_Al_RelPos),]
dfFilt$Var_Al_RelPos = as.numeric(as.character(dfFilt$Var_Al_RelPos))
dfFilt$Ref_Al_RelPos = as.numeric(as.character(dfFilt$Ref_Al_RelPos))

model1 <- glm(ConsTest ~ ALLELE.FREQUENCY+STRAND.BIAS+DEPTH+QUAL+Var_Al_RelPos+Ref_Al_RelPos, data = dfFilt, family = "binomial")
summary(model1)

step_model  = stepAIC(model1, direction = "both" , trace = T, steps = 10000)
summary(step_model) # ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos

sumModel = function(model) {
  modSum = summary(model)
  modRes = data.frame(modSum[["coefficients"]])
  colnames(modRes)[4] = "Pval"
  modRes$Variables = rownames(modRes)
  modRes$AIC = modSum$aic
  return(modRes)
}

multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos, data = dfFilt, family = "binomial")

modRes2 = sumModel(multivarModel)

write.csv(modRes2, "metaseq_ViralrecLofreq_coeficients.csv", row.names = F)

system("aws s3 cp metaseq_ViralrecLofreq_coeficients.csv s3://abombin/Vivacity/classifier/")

#AUC
set.seed(123)
train_idx <- createDataPartition(dfFilt$ConsTest, p = 0.7, list = FALSE)
train_data <- dfFilt[train_idx, ]
test_data <- dfFilt[-train_idx, ]
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos, data = dfFilt, family = "binomial")
probs <- predict(multivarModel, newdata = test_data, type = "response")
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)