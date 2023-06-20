import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

df = pd.read_csv('Ludy_metaAmpIvar_overlapSnv_RelPos.csv')
mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

len(varNa.index)

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)

# best on auc relPos_2 has the highest score 99.65
best_aic2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV']
relPos = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV', 'Var_Al_RelPos', 'Ref_Al_RelPos']
best_aic3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos']
best_aic3_mod1 = ['ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos']
best_aic3_mod2 = ['ALT_FREQ', 'ALT_DP', 'REF_DP', 'ALT_RV', 'Var_Al_RelPos']

X = dfFilt[best_aic3_mod1]
y = dfFilt['ConsTest']

#varNa = X[X.isna().any(axis=1)]
#respNa = y[y.isna()]

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Create and train the logistic regression model
model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")
