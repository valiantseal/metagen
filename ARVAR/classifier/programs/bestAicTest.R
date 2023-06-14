library(MASS)
df = read.csv("288_metaAmpIvar_overlapSnv.csv")
df = read.csv("Ludy_metaAmpIvar_overlapSnv.csv")

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

# previous aic model would be wrong since several variables will be clearly linearly related
# better one would be ConsTest ~ ALT_FREQ + ALT_QUAL + ALT_DP + REF_DP + REF_QUAL + ALT_RV
# for 288 dataset best model is ALT_FREQ + ALT_QUAL + ALT_DP + REF_RV + ALT_RV