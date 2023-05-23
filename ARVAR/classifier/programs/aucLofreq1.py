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
    

sumStat = Summary(dfFilt = dfMetaFilt, perSample = True)
sumIvar = Summary(dfFilt = dfRef, perSample = True)
sumIvarSub = sumIvar[["Sample", "True_Positive"]]
sumIvarSub.columns = ["Sample", "Total_Positive_Truth"]
combSummary = pd.merge(sumStat, sumIvarSub, how = "left", on = "Sample")
combSummary.to_csv(path_or_buf = "test_consensus/LofreqIvarSummary.csv", index = False)


# Assuming you have a pandas DataFrame with your data
# X represents the features, and y represents the target variable
colOpt1 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos', 'Ref_Al_RelPos']
colOpt2 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos']
colOpt5 = ['ALLELE.FREQUENCY', 'STRAND.BIAS', 'DEPTH', 'QUAL', 'Var_Al_RelPos','Ref_Al_RelPos', "Sample", 'Samp_Pos_Ref_Alt']

X = dfMetaFilt[colOpt5]
y = dfMetaFilt['ConsTest']
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
  dfFilt = df[(df["ALLELE.FREQUENCY"] >= minFreq) & (df["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
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
  dfFilt = df[(df["ALLELE.FREQUENCY"] >= minFreq) & (df["ALLELE.FREQUENCY"] <= maxFreq)].reset_index().drop(["index"], axis =1)
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

sumTest.to_csv(path_or_buf = "test_consensus/LofreqIvarTestSummary.csv", index = False)
