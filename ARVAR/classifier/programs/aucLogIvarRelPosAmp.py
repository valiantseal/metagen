import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score


df = pd.read_csv('Ludy_ampDedupMetaIvar_overlapSnv_RelPos.csv')
mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)

best_aic = ["ALT_FREQ", "ALT_QUAL", "REF_QUAL", "Ref_Al_RelPos"] # 0.9148099318831027
best_aic3_mod1 = ['ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos'] # 0.9048209184794551
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV'] # 0.9125994286969896  fullDf = 0.9294278133760512
allVars = ['ALT_FREQ','ALT_QUAL','ALT_DP','REF_DP','REF_QUAL','REF_RV','ALT_RV','Var_Al_RelPos','Ref_Al_RelPos'] # 0.91121072291804
best_aic_orig = ["ALT_FREQ" , "ALT_QUAL", "REF_QUAL", "REF_RV"] # 0.9114590199956055  fullDf = 0.9288229208383393
best_aic_var = ["ALT_FREQ", "ALT_QUAL", "REF_QUAL", "Ref_Al_RelPos", "Var_Al_RelPos"] # 0.9115776752362118

#X = dfFilt[best_aic_var]
#y = dfFilt['ConsTest']

X = df[best_aic_orig]
y = df['ConsTest']

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
