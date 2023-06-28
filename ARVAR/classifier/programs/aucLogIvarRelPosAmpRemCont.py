import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


df = pd.read_csv('result_tables/Ludy_ampDedupMetaIvar_overlapSnv_RelPos.csv')
mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)
exclList = ['EHC-C19-2346H', 'EHC-C19-2387W', 'EHC-C19-2777U']
exclDf = dfFilt[dfFilt["Sample"].isin(exclList)]
dfFilt = dfFilt[~dfFilt["Sample"].isin(exclList)]

dfFilt.to_csv("result_tables/Ludy_ampDedupMetaIvar_overlapSnv_RelPosRemCont.csv", index = False)

best_aic = ["ALT_FREQ", "ALT_QUAL", "REF_QUAL", "Ref_Al_RelPos"] # 0.9191089883295839
best_aic1 = ["ALT_FREQ", "ALT_QUAL", "REF_QUAL", "REF_RV", "Var_Al_RelPos"] # 0.9043012689820022
best_aic2 = ['ALT_FREQ', 'ALT_QUAL', 'REF_QUAL', 'REF_RV'] # 0.9146271266872892

X = dfFilt[best_aic2]
y = dfFilt['ConsTest']



# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Create and train the logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

coefficients = model.coef_[0]

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")
