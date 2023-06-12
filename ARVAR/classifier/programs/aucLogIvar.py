from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALT_FREQ"] >= minFreq) & (metaseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  ampseqFilt =  ampseq[(ampseq["ALT_FREQ"] >= minFreq) & (ampseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
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
dfFilt = getConsensus(metaSeq = "test_consensus/metaseqIvar.csv", ampSeq = "test_consensus/ampseqConsIvar.csv", minFreq = 0.02, maxFreq = 1)

#dfFilt.to_csv("Ludy_metaAmpIvar_overlapSnv.csv", index = False)

value_counts = dfFilt['ConsTest'].value_counts()
print(value_counts)

def Summary(dfFilt, perSample):
  snvs = []
  truePos = []
  falsePos = []
  Samples = []
  if perSample == True:
    samplesList = list(pd.unique(dfFilt["Sample"]))
    for sample in samplesList:
      dfSub = dfFilt[dfFilt["Sample"] == sample]
      dfPos = dfSub[dfSub["ConsTest"] ==1].reset_index().drop(["index"], axis =1)
      dfNeg = dfSub[dfSub["ConsTest"] ==0].reset_index().drop(["index"], axis =1)
      curSnvs = dfSub['Samp_Pos_Ref_Alt'].nunique()
      curPos = dfPos['Samp_Pos_Ref_Alt'].nunique()
      curNeg = dfNeg['Samp_Pos_Ref_Alt'].nunique()
      
      snvs.append(curSnvs)
      truePos.append(curPos)
      falsePos.append(curNeg)
      Samples.append(sample)
    combDat = pd.DataFrame({"Sample": Samples, "Total_SNVs": snvs, "True_Positive": truePos, 'False_Positive': falsePos})
    return(combDat)
      
  elif perSample == False:
    dfPos = dfFilt[dfFilt["ConsTest"] ==1].reset_index().drop(["index"], axis =1)
    dfNeg = dfFilt[dfFilt["ConsTest"] ==0].reset_index().drop(["index"], axis =1)
    print('Total number of unique SNVs  ' + str(dfFilt['Samp_Pos_Ref_Alt'].nunique()))
    print('Number of true positive SNVs  ' + str(dfPos['Samp_Pos_Ref_Alt'].nunique()))
    print('Number of false positive SNVs  ' + str(dfNeg ['Samp_Pos_Ref_Alt'].nunique()))
    

sumStat = Summary(dfFilt = dfFilt, perSample = True)

#sumStat.to_csv(path_or_buf = "test_consensus/IvarSummary.csv", index = False)

list(dfFilt)
# Assuming you have a pandas DataFrame with your data
# X represents the features, and y represents the target variable
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']
colOpt2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL']
colOpt3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP']
colOpt4 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'ALT_RV']
colOpt5 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'Samp_Pos_Ref_Alt', "Sample"]

X = dfFilt[colOpt1]
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

## alternative version to extract predictions
X = dfFilt[colOpt5]
y = dfFilt['ConsTest']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

X_train1 = X_train[colOpt1]
X_test1 = X_test[colOpt1]

model = LogisticRegression(max_iter=1000)
model.fit(X_train1, y_train)
y_pred_proba = model.predict_proba(X_test1)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")


val_predict=model.predict(X_test1)
# make data with Ids and survival
X_test["Predict_val"] = val_predict
X_test['ConsTest'] = y_test

def SumPredict(df, minFreq, maxFreq):
  corPredict = []
  wrongPredict = []
  dfFilt = df[(df["ALT_FREQ"] >= minFreq) & (df["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  for i in range(len(dfFilt.index)):
    if dfFilt.loc[i, "Predict_val"] == dfFilt.loc[i, "ConsTest"]:
      corPredict.append(dfFilt.loc[i, "Samp_Pos_Ref_Alt"])
    else:
      wrongPredict.append(dfFilt.loc[i, "Samp_Pos_Ref_Alt"])
  numbCorPred = len(corPredict)
  numbWrongPred = len(wrongPredict)
  return numbCorPred, numbWrongPred

numbCorPred, numbWrongPred = SumPredict(df = X_test, minFreq = 0.02, maxFreq = 1)

numbCorPred / (numbCorPred + numbWrongPred) * 100
  
def SumPredictPerSamp(df, minFreq, maxFreq):
  snvs = []
  truePos = []
  falsePos = []
  Samples = []
  totalPos = []
  corIdent = []
  wrongIdnet = []
  falseNeg = []
  dfFilt = df[(df["ALT_FREQ"] >= minFreq) & (df["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  samples = list(pd.unique(df["Sample"]))
  for sample in samples:
    corPredict = []
    wrongPredict = []
    curTruePos = []
    curFalsePos = []
    curFalseNeg = []
    dfSub = dfFilt[dfFilt["Sample"] == sample].reset_index().drop(["index"], axis =1)
    for i in range(len(dfSub.index)):
      if dfSub.loc[i, "Predict_val"] == dfSub.loc[i, "ConsTest"]:
        corPredict.append(dfSub.loc[i, "Samp_Pos_Ref_Alt"])
      else:
        wrongPredict.append(dfSub.loc[i, "Samp_Pos_Ref_Alt"])
        if dfSub.loc[i, "Predict_val"] == 1:
          curFalsePos.append(dfSub.loc[i, "Samp_Pos_Ref_Alt"])
        elif dfSub.loc[i, "Predict_val"] == 0:
          curFalseNeg.append(dfSub.loc[i, "Samp_Pos_Ref_Alt"])
    numbCorPred = len(corPredict)
    numbWrongPred = len(wrongPredict)
    dfPos = dfSub[dfSub["ConsTest"] ==1].reset_index().drop(["index"], axis =1)
    curPos = dfPos['Samp_Pos_Ref_Alt'].nunique()
    curSnvs = dfSub['Samp_Pos_Ref_Alt'].nunique()
    for i in range(len(dfPos.index)):
      if dfPos.loc[i, "Predict_val"] == dfPos.loc[i, "ConsTest"]:
        curTruePos.append(dfPos.loc[i, "Samp_Pos_Ref_Alt"])
        
    numbTruePos = len(curTruePos)
    numbFalsePos = len(curFalsePos)
    numbFalseNeg = len(curFalseNeg)
    
    #append values
    snvs.append(curSnvs)
    corIdent.append(numbCorPred)
    wrongIdnet.append(numbWrongPred)
    Samples.append(sample)
    totalPos.append(curPos)
    truePos.append(numbTruePos)
    falsePos.append(numbFalsePos)
    falseNeg.append(numbFalseNeg)
  combDat = pd.DataFrame({"Sample": Samples, "Total_SNVs": snvs, "Correctly_identified": corIdent, "Incorrectly_identified": wrongIdnet ,
  "True_Positive": truePos, 'False_Positive': falsePos, "False_Negative":falseNeg, "Total_Positive_Truth": totalPos})
  return combDat

sumTest = SumPredictPerSamp(df = X_test, minFreq = 0.02, maxFreq = 1)

#sumTest.to_csv(path_or_buf = "test_consensus/IvarTestSummary.csv", index = False)
