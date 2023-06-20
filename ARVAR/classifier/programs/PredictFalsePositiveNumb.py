X_train = dfFilt[colOpt1].reset_index(drop = True)
y_train = dfFilt['ConsTest']


X_test = X_train
y_test = y_train

X_test1 = dfFilt[colOpt5].reset_index(drop = True)

model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)
y_pred_proba = model.predict_proba(X_train)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_train, y_pred_proba)

print(f"AUC Score: {auc_score}")

val_predict=model.predict(X_test)

# make data with Ids and survival
X_test1["Predict_val"] = val_predict
X_test1['ConsTest'] = y_test

consSnv = X_test1[X_test1["ConsTest"] == 1]

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

sumTest = SumPredictPerSamp(df = X_test1, minFreq = 0.02, maxFreq = 1)

len(pd.unique(X_test1["Samp_Pos_Ref_Alt"]))

column_sum = sumTest['False_Positive'].sum()
column_std = sumTest['False_Positive'].std()

column_sum/len(pd.unique(X_test1["Samp_Pos_Ref_Alt"])) 

0.00399814039981404 * 6220

###
len(pd.unique(X_test["Samp_Pos_Ref_Alt"]))

column_sum = sumTest['False_Positive'].sum()
column_std = sumTest['False_Positive'].std()

column_sum/len(pd.unique(X_test["Samp_Pos_Ref_Alt"])) 

0.0037186241090796405 * 6220
