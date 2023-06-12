df = read.csv("Ludy_metaAmpIvar_overlapSnv.csv")


model1 <- glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_QUAL + ALT_RV + TOTAL_DP, data = df, family = "binomial")

summary(model1)

model2 =  glm(ConsTest ~ ALT_FREQ + REF_DP + REF_QUAL + ALT_RV, data = df, family = "binomial")

summary(model2)


model2 =  glm(ConsTest ~ ALT_QUAL, data = df, family = "binomial")

summary(model2)

# get AUC
library(pROC)
library(caret)


set.seed(123)
train_idx <- createDataPartition(df$ConsTest, p = 0.7, list = FALSE)
train_data <- df[train_idx, ]
test_data <- df[-train_idx, ]

multivarModel =  glm(ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_QUAL + ALT_RV + TOTAL_DP, data = df, family = "binomial")
summary(multivarModel)

probs <- predict(multivarModel, newdata = test_data, type = "response")

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)