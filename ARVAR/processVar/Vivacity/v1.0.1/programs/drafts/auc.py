from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd

df = pd.read_csv("test_consensus/vivacityNoFilt_ivarfreq0.02.csv")

def filterDat(df, minFreq):
  dfFilt = df.dropna(subset=['Var_Al_RelPos'])
  dfFreq = dfFilt[dfFilt["ALLELE.FREQUENCY"] >= minFreq]
  return dfFreq

dfFilt = filterDat(df = df, minFreq = 0.01)

# Assuming you have a pandas DataFrame with your data
# X represents the features, and y represents the target variable
X = dfFilt[['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos']]
y = dfFilt['ConsTest']

q1 = dfFilt[dfFilt.isna().any(axis=1)]

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=13)

# Create and train the logistic regression model
model = LogisticRegression()
model.fit(X_train, y_train)

# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")
