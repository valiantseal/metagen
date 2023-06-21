from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os
from itertools import combinations

allVars = ['ALT_FREQ','ALT_QUAL','ALT_DP','REF_DP','REF_QUAL','REF_RV','ALT_RV','Var_Al_RelPos','Ref_Al_RelPos']

def getBestAuc(df , allVars):
  best_auc_score = 0.0
  best_features = []
  combFet = []
  combAuc = []
  curDf = pd.read_csv(df)
  dfFilt = curDf.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)
  X = dfFilt[allVars]
  y = dfFilt['ConsTest']
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
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
        
        # Update the best AUC score and feature combination if necessary
        if auc_score > best_auc_score:
            best_auc_score = auc_score
            best_features = selected_features
            
  combDat = pd.DataFrame({"Features":combFet, "AUC":combAuc})

  # Print the best feature combination and AUC score
  print("Best feature combination:", best_features)
  print("Best AUC score:", best_auc_score)
  return  combDat

aucMeta = getBestAuc(df = 'Ludy_metaAmpIvar_overlapSnv_RelPos.csv' , allVars = allVars) # 'ALT_FREQ', 'ALT_QUAL', 'REF_QUAL', 'REF_RV', 'Var_Al_RelPos'
aucAmp = getBestAuc(df = "Ludy_ampDedupMetaIvar_overlapSnv_RelPos.csv" , allVars = allVars) # 'ALT_FREQ', 'ALT_DP', 'REF_QUAL', 'ALT_RV', 'Ref_Al_RelPos'
       
