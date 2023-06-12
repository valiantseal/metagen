import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os


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

best_auc = ['ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV']
best_aic = ['ALT_FREQ' , 'ALT_QUAL' , 'ALT_DP' , 'REF_QUAL' , 'ALT_RV' , 'TOTAL_DP']
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']

X = dfFilt[best_auc]
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

coef_df.to_csv("LogIvarMeta_bestAuc_CoefPval.csv", index = False)

# best AIC
X = dfFilt[best_aic]
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

coef_df.to_csv("LogIvarMeta_bestAic_CoefPval.csv", index = False)

# model that was currently used
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

coef_df.to_csv("LogIvarMeta_CurUsed_CoefPval.csv", index = False)
