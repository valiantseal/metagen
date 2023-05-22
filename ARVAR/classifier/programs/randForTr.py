from sklearn.ensemble import RandomForestClassifier

colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos']
colOpt2 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos']

X = dfMetaFilt[colOpt1]
y = dfMetaFilt['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

model = RandomForestClassifier(criterion='gini', 
                             n_estimators=1200,
                             min_samples_split=10,
                             min_samples_leaf=1,
                             max_features='auto',
                             oob_score=True,
                             random_state=42,
                             n_jobs=-1)
model.fit(X_train, y_train)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")


print("%.4f" % model.oob_score_)
# show variables importance 
pd.concat((pd.DataFrame(X_train.columns, columns = ['variable']), 
           pd.DataFrame(model.feature_importances_, columns = ['importance'])), 
          axis = 1).sort_values(by='importance', ascending = False)[:20]
