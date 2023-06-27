from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALLELE.FREQUENCY"] >= minFreq) & (metaseq["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  ampseqFilt =  ampseq[(ampseq["ALLELE.FREQUENCY"] >= minFreq) & (ampseq["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  targSnv = ampseqFilt["Samp_Pos_Ref_Alt"].to_list()
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
df = getConsensus(metaSeq = "test_consensus/metaseq_comb_vivacity.csv", ampSeq = "test_consensus/ampseq_comb_vivacity.csv", minFreq = 0.02, maxFreq = 1)
dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)

exclList = ['EHC-C19-2346H-L2', 'EHC-C19-2387W-L2', 'EHC-C19-2777U-L2']
dfFilt = dfFilt[~dfFilt["ExactSample"].isin(exclList)]

dfFilt.to_csv("result_tables/Ludy_metaAmp_ViralrecLofreq_overlapSnv_RemCont.csv", index = False)

value_counts = dfFilt['ConsTest'].value_counts()
print(value_counts)

colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos'] # 0.9297619867336868
colOpt5 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos','Ref_Al_RelPos', "Sample", 'Samp_Pos_Ref_Alt']
best_aic = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos'] # 0.948325608294075

X = dfFilt[best_aic]
y = dfFilt['ConsTest']

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
