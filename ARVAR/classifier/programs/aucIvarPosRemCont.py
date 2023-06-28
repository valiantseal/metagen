import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

df = pd.read_csv('result_tables/Ludy_metaAmpIvar_overlapSnv_RelPos.csv')
mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

len(varNa.index)

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)
'EHC-C19-2777U' in dfFilt["Sample"].to_list()
'EHC-C19-2346H' in dfFilt["Sample"].to_list()
'EHC-C19-2387W' in dfFilt["Sample"].to_list()
exclList = ['EHC-C19-2346H-L2', 'EHC-C19-2387W-L2', 'EHC-C19-2777U-L2']

exclDf = dfFilt[dfFilt["ExactSamp"].isin(exclList)]

dfFilt = dfFilt[~dfFilt["ExactSamp"].isin(exclList)]
'EHC-C19-2985U' in dfFilt["Sample"].to_list()

blackList = pd.read_csv("samples.blacklist", header=None, names=["ExactSamp"])
blackList = blackList["ExactSamp"].to_list()
black_list = []
for i in blackList:
  newSamp = i.replace("_", "-")
  black_list.append(newSamp)

blackDf = dfFilt[dfFilt["ExactSamp"].isin(black_list)]
pd.unique(blackDf["ExactSamp"])

dfFilt = dfFilt[~dfFilt["ExactSamp"].isin(black_list)]
pd.unique(dfFilt["ExactSamp"])

#dfFilt.to_csv("Ludy_metaAmpIvar_overlapSnv_RelPos_RemCont.csv", index = False)

# best on auc relPos_2 has the highest score 99.65
best_aic2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV']
relPos = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV', 'Var_Al_RelPos', 'Ref_Al_RelPos']
best_aic3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos'] # 0.9983879482775286
best_aic3_mod1 = ['ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos']
best_aic3_mod2 = ['ALT_FREQ', 'ALT_DP', 'REF_DP', 'ALT_RV', 'Var_Al_RelPos']
best_aic_remCont = ["ALT_FREQ", "ALT_DP", "REF_DP", "REF_QUAL", "ALT_RV", "Var_Al_RelPos"] # 0.9984444593157702

X = dfFilt[best_aic_remCont]
y = dfFilt['ConsTest']

#varNa = X[X.isna().any(axis=1)]
#respNa = y[y.isna()]

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Create and train the logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

coefficients = model.coef_[0]
feature_names = best_aic_remCont

# Pair the coefficients with feature names
coef_with_names = list(zip(feature_names, coefficients))

# Print the coefficients with their names
for name, coef in coef_with_names:
    print(name, coef)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")
