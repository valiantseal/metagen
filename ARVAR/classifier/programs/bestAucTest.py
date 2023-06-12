from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
from itertools import combinations

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

# Assuming you have a pandas DataFrame with your data
# X represents the features, and y represents the target variable
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']
colOpt2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL']
colOpt3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP']
colOpt4 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'ALT_RV']
colOpt5 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'Samp_Pos_Ref_Alt', "Sample"]
colOpt6 = ['ALT_FREQ' , 'ALT_QUAL' , 'ALT_DP' , 'REF_QUAL' , 'ALT_RV' , 'TOTAL_DP'] # best AIC model
X = dfFilt[colOpt3]
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


####
best_auc_score = 0.0
best_features = []

combFet = []
combAuc = []

# Loop through all possible combinations of feature indices
for r in range(1, len(X.columns) + 1):
    for feature_indices in combinations(range(len(X.columns)), r):
        # Select the corresponding features from the DataFrame
        selected_features = X.columns[list(feature_indices)]
        
        # Train a logistic regression model using the selected features
        model = LogisticRegression(max_iter=1000)
        model.fit(X_train[selected_features], y_train)
        
        # Predict probabilities for the test data
        y_pred_proba = model.predict_proba(X_test[selected_features])[:, 1]
        
        # Calculate the AUC score
        auc_score = roc_auc_score(y_test, y_pred_proba)
        
        curFeature = '--'.join(selected_features)
        curAuc = auc_score
        
        combFet.append(curFeature)
        combAuc.append(curAuc)
        
        combDat = pd.DataFrame({"Features":combFet, "AUC":combAuc})
        
        # Update the best AUC score and feature combination if necessary
        if auc_score > best_auc_score:
            best_auc_score = auc_score
            best_features = selected_features

# Print the best feature combination and AUC score
print("Best feature combination:", best_features)
print("Best AUC score:", best_auc_score)

# beas AUC model ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV'
# ALT_FREQ--REF_DP--REF_QUAL--ALT_RV

# if use on all data without splitting in training and testing 'ALT_FREQ', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'ALT_RV'
