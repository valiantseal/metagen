from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd


minFreq = 0.02
maxFreq = 1

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALT_FREQ"] >= minFreq) & (metaseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  ampseqFilt =  ampseq[(ampseq["ALT_FREQ"] >= minFreq) & (ampseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  targSnv = ampseq["Samp_Pos_Ref_Alt"].to_list()
  ConsTest = []
  for i in range(len(metaseqFilt.index)):
    curSnv = metaseqFilt.loc[i, "Samp_Pos_Ref_Alt"]
    if curSnv in targSnv:
      curCons = 1
    else:
      curCons = 0
    ConsTest.append(curCons)
  metaseqFilt["ConsTest"] = ConsTest
  return(metaseqFilt)

# increasing lower frequency increase accuracy 
# colOption1, minFreq at least 0.02 and ampSeqIvar have the bset performance   
dfFilt = getConsensus(metaSeq = "test_consensus/metaseqIvar.csv", ampSeq = "test_consensus/ampseqConsIvar.csv", minFreq = minFreq, maxFreq = maxFreq)
dfRef = dfFilt[dfFilt["ConsTest"] ==1].reset_index().drop(["index"], axis =1)

# get metaseq

dfMeta = pd.read_csv("test_consensus/metaseqLofreq.csv")

def annotLofreq(df, ref, minFreq, maxFreq):
  ConsTest = []
  metaseqFilt = df[(df["ALLELE.FREQUENCY"] >= minFreq) & (df["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  metaseqFilt = metaseqFilt.dropna(subset=['Var_Al_RelPos']).reset_index().drop(["index"], axis =1)
  targSnv = ref["Samp_Pos_Ref_Alt"].to_list()
  for i in range(len(metaseqFilt.index)):
    curSnv = metaseqFilt.loc[i, "Samp_Pos_Ref_Alt"]
    if curSnv in targSnv:
      curCons = 1
    else:
      curCons = 0
    ConsTest.append(curCons)
  metaseqFilt["ConsTest"] = ConsTest
  return(metaseqFilt)

dfMetaFilt = annotLofreq(df = dfMeta, ref = dfRef, minFreq = minFreq, maxFreq = maxFreq)
value_counts = dfMetaFilt['ConsTest'].value_counts()
print(value_counts)

q1 = dfMetaFilt[dfMetaFilt.isna().any(axis=1)]

list(dfMetaFilt)
# Assuming you have a pandas DataFrame with your data
# X represents the features, and y represents the target variable
colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos']
colOpt2 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos']

X = dfMetaFilt[colOpt2]
y = dfMetaFilt['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

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

