import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
import subprocess

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

dfFilt = getConsensus(metaSeq = "test_consensus/288_metaseqIvar.csv", ampSeq = "test_consensus/288_ampseqIvar.csv", minFreq = 0.02, maxFreq = 1)

best_aic2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV']
best_aic_288 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_RV', 'ALT_RV']
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']

X = dfFilt[best_aic2]
y = dfFilt['ConsTest']
logit_model=sm.Logit(y,X)
result=logit_model.fit()
print(result.summary())
summary_table = result.summary2()
# Extract the coefficients and p-values
coef_df = summary_table.tables[1]
coef_df["Variable"] = coef_df.index
# Print the coefficient and p-value DataFrame
print(coef_df)
coef_df.to_csv("288Meta_bestAic2_CoefPval.csv", index = False)

# best aic specific to 288 samples
X = dfFilt[best_aic_288]
y = dfFilt['ConsTest']
logit_model=sm.Logit(y,X)
result=logit_model.fit()
print(result.summary())
summary_table = result.summary2()
# Extract the coefficients and p-values
coef_df = summary_table.tables[1]
coef_df["Variable"] = coef_df.index
# Print the coefficient and p-value DataFrame
print(coef_df)
coef_df.to_csv("288Meta_best288Aic_CoefPval.csv", index = False)

# currently used model
X = dfFilt[colOpt1]
y = dfFilt['ConsTest']
logit_model=sm.Logit(y,X)
result=logit_model.fit()
print(result.summary())
summary_table = result.summary2()
# Extract the coefficients and p-values
coef_df = summary_table.tables[1]
coef_df["Variable"] = coef_df.index
# Print the coefficient and p-value DataFrame
print(coef_df)
coef_df.to_csv("288Meta_CurUsed_CoefPval.csv", index = False)

# transfer to S3
def toS3(filesList):
  for i in filesList:
    try:
      cmd_str = f'aws s3 cp {i} s3://abombin/Vivacity/classifier/'
      subprocess.run(cmd_str, shell = True)
    except:
      pass
    
filesList = ["288Meta_bestAic2_CoefPval.csv", "288Meta_best288Aic_CoefPval.csv", "288Meta_CurUsed_CoefPval.csv"]

toS3(filesList = filesList)
