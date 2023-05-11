library(caret)
library(pROC)


# Split data into training and test sets
set.seed(123)
train_idx <- createDataPartition(metaResDf$ConsTest, p = 0.7, list = FALSE)
train_data <- metaResDf[train_idx, ]
test_data <- metaResDf[-train_idx, ]

# Fit a logistic regression model on the training set
model <- glm(ConsTest ~ Position_test, 
             data = train_data, family = binomial)

# Get predicted probabilities for the test set
probs <- predict(model, newdata = test_data, type = "response")

# Determine optimal threshold using AUC
roc_obj <- roc(test_data$ConsTest, probs)
roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

plot(roc_obj)

multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + DEPTH + QUAL + Position_test, data = train_data, family = binomial)

probs <- predict(multivarModel, newdata = test_data, type = "response")

roc_obj <- roc(test_data$ConsTest, probs)

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)

plot(roc_obj)

test_data$Prob = predict(multivarModel, newdata = test_data, type = "response")

# most reasonable model without exact value for position test
# Split data into training and test sets
set.seed(123)
train_idx <- createDataPartition(metaResDf$ConsTest, p = 0.7, list = FALSE)
train_data <- metaResDf[train_idx, ]
test_data <- metaResDf[-train_idx, ]

multivarModel = glm(ConsTest ~ RawVarFreq + STRAND.BIAS + DEPTH + QUAL + Var_Al_RelPos, data = train_data, family = binomial)
summary(multivarModel)

probs <- predict(multivarModel, newdata = test_data, type = "response")

roc_obj <- roc(test_data$ConsTest ~ probs, plot = TRUE, print.auc = TRUE)