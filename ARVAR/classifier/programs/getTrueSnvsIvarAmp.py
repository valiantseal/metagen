from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALT_FREQ"] >= minFreq) & (metaseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  filter_list = metaseqFilt["Sample"].to_list()
  ampseqFilt =  ampseq[(ampseq["ALT_FREQ"] >= minFreq) & (ampseq["ALT_FREQ"] <= maxFreq) & (ampseq['Sample'].isin(filter_list))].reset_index().drop(["index"], axis =1)
  targSnv = metaseqFilt["Samp_Pos_Ref_Alt"].to_list()
  ConsTest = []
  for i in range(len(ampseqFilt.index)):
    curSnv = ampseqFilt.loc[i, "Samp_Pos_Ref_Alt"]
    if curSnv in targSnv:
      curCons = 1
    else:
      curCons = 0
    ConsTest.append(curCons)
  ampseqFilt["ConsTest"] = ConsTest
  return(ampseqFilt)

# increasing lower frequency increase accuracy 
# colOption1, minFreq at least 0.02 and ampSeqIvar have the bset performance   
dfFilt = getConsensus(metaSeq = "test_consensus/metaseqIvar.csv", ampSeq = "test_consensus/ampseqIvar.csv", minFreq = 0.02, maxFreq = 1)

colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']
colOpt2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL']
colOpt3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP']
colOpt4 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'ALT_RV']
colOpt5 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'Samp_Pos_Ref_Alt', "Sample"]

X_train = dfFilt[colOpt1]
y_train = dfFilt['ConsTest']


model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)
# Predict probabilities on the test set
y_pred_proba = model.predict_proba(X_train)[:, 1]
# Calculate the AUC score
auc_score = roc_auc_score(y_train, y_pred_proba)
print(f"AUC Score: {auc_score}")

curSamples = list(pd.unique(dfFilt["Sample"]))

# get Ludy's data
def combIvarDat(path, minFreq, maxFreq, curSamples):
  filesList = os.listdir(path)
  combDat = pd.DataFrame()
  for i in filesList:
    inFile = path + i
    try:
      curDf = pd.read_table(inFile, low_memory=False)
      dfFilt = curDf[(curDf['PASS'] == True) & (curDf['ALT_FREQ'] >= minFreq) & (curDf['ALT_FREQ'] <= maxFreq)].reset_index().drop(["index"], axis =1)
      if len(dfFilt.index) < 2:
        print(inFile + "no  SNVs")
      sampleName = i.replace(".ivar_trim.sorted.bam_snps.tsv", "")
      dfFilt["ExactSamp"] = sampleName
      sampleList = sampleName.split("_")[0]
      sampleList = sampleList.split("-")[0:3]
      newName = "-".join(sampleList)
      dfFilt["Sample"] = newName 
      filtered_df = dfFilt[~dfFilt['Sample'].isin(curSamples)].reset_index().drop(["index"], axis =1)
      combDat = pd.concat([combDat, filtered_df])
    except Exception as e:
      print(i)
      print(e)
  combDat = combDat.reset_index().drop(["index"], axis =1)
  return combDat

ivarDat = combIvarDat(path = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/relaxIvar/",
                      minFreq = 0.02, maxFreq = 1, curSamples = curSamples)
                      
len(pd.unique(ivarDat["ExactSamp"]))
len(pd.unique(dfFilt["FullSamp"]))

def addUnSnv(df):
  Samp_Pos_Ref_Alt= []
  for i in (range(len(df.index))):
    curList = [df.loc[i, 'Sample'], str(df.loc[i, 'POS']), df.loc[i, 'REF'], df.loc[i, 'ALT']]
    curSamp = "__".join(curList)
    Samp_Pos_Ref_Alt.append(curSamp)
  df["Samp_Pos_Ref_Alt"] = Samp_Pos_Ref_Alt
  return df

ivarDat = addUnSnv(df = ivarDat)    

# predict
X_test = ivarDat[colOpt1]

val_predict=model.predict(X_test)
# make data with Ids and survival
ivarDat["ConsTest"] = val_predict

def combTrueSnv(df1, df2):
  df1Filt = df1[df1["ConsTest"] == 1].reset_index().drop(["index"], axis =1)
  df2Filt = df2[df2["ConsTest"] == 1].reset_index().drop(["index"], axis =1)
  combDat = pd.concat([df1Filt, df2Filt]).reset_index().drop(["index"], axis =1)
  return combDat

trueSnvs = combTrueSnv(df1 = dfFilt, df2 = ivarDat)

trueSnvs.to_csv("Ludy_AmpSeq_TrueSnvs_LogClass.csv", index = False)
