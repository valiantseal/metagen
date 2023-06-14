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

colOpt3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP']
colOpt4 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'ALT_RV']
colOpt5 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'Samp_Pos_Ref_Alt', "Sample"]
colOpt6 = ['ALT_FREQ' , 'ALT_QUAL' , 'ALT_DP' , 'REF_QUAL' , 'ALT_RV' , 'TOTAL_DP'] # best AIC model
X = dfFilt[colOpt3]
y = dfFilt['ConsTest']

combFet = []
combAuc = []

# Loop through all possible combinations of feature indices
for r in range(1, len(X.columns) + 1):
    for feature_indices in combinations(range(len(X.columns)), r):
        # Select the corresponding features from the DataFrame
        selected_features = X.columns[list(feature_indices)]
        
        # Train a logistic regression model using the selected features
        model = LogisticRegression(max_iter=1000)
        model.fit(X[selected_features], y)
        
        logit_model=sm.Logit(y,X[selected_features])
        result=logit_model.fit()

        summary_table = result.summary2()

        # Extract the coefficients and p-values
        coef_df = summary_table.tables[0]
        aic = coef_df.iloc[2,3]
        
        
        curFeature = '--'.join(selected_features)
        curAic = aic
        
        combFet.append(curFeature)
        combAuc.append(aic)
        
        combDat = pd.DataFrame({"Features":combFet, "AIC":combAuc})
        
combDat['AIC'] = combDat['AIC'].astype(float)
