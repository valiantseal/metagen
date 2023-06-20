import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
import pandas as pd
import os

df = pd.read_csv('Ludy_metaAmpIvar_overlapSnv_RelPos.csv')
mask = df[['Var_Al_RelPos', 'Ref_Al_RelPos']].isna().any(axis=1)
varNa = df[mask]

len(varNa.index)

dfFilt = df.dropna(subset=['Var_Al_RelPos', 'Ref_Al_RelPos'], how='any').reset_index(drop = True)

# best on auc relPos_2 has the highest score 99.65
best_aic2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV']
relPos = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP' , 'REF_QUAL' , 'ALT_RV', 'Var_Al_RelPos', 'Ref_Al_RelPos']
best_aic3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos']
best_aic3_mod1 = ['ALT_FREQ', 'REF_DP', 'REF_QUAL', 'ALT_RV', 'Var_Al_RelPos']
best_aic3_mod2 = ['ALT_FREQ', 'ALT_DP', 'REF_DP', 'ALT_RV', 'Var_Al_RelPos']

# best AIC2
X = dfFilt[best_aic3_mod2]
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
#coef_df.to_csv("LogIvarMeta_BestAic2_CoefPval.csv", index = False)
