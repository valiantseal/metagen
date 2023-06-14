from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
import pandas as pd
import matplotlib.pyplot as plt

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
#dfFilt.to_csv("288_metaAmpIvar_overlapSnv.csv", index = False)

value_counts = dfFilt['ConsTest'].value_counts()
print(value_counts)

# X represents the features, and y represents the target variable
colOpt1 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV']
colOpt2 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL']
colOpt3 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'REF_DP', 'REF_QUAL', 'REF_RV', 'ALT_RV', 'TOTAL_DP']
colOpt4 = ['ALT_FREQ', 'ALT_QUAL', 'ALT_DP', 'ALT_RV']


X = dfFilt[colOpt1]
y = dfFilt['ConsTest']

varNa = X[X.isna().any(axis=1)]
respNa = y[y.isna()]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)


model = LogisticRegression(max_iter=1000)
model.fit(X_train, y_train)
y_pred_proba = model.predict_proba(X_test)[:, 1]

# Calculate the AUC score
auc_score = roc_auc_score(y_test, y_pred_proba)

print(f"AUC Score: {auc_score}")

# plot

def plot_auc_curve(y_test, y_pred_proba):
    # Create a new figure and axes
    fig, ax = plt.subplots(figsize=(10, 8))  # Adjust width and height as desired

    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)

    # Clear the previous plot
    ax.clear()

    # Plotting the AUC curve
    ax.plot(fpr, tpr, label='AUC Curve (AUC = %0.2f)' % auc_score)
    ax.plot([0, 1], [0, 1], 'k--', label='Random')
    ax.set_xlabel('False Positive Rate (FPR)')
    ax.set_ylabel('True Positive Rate (TPR)')
    ax.set_title('Receiver Operating Characteristic (ROC)')
    ax.legend(loc='lower right')

    # Save the plot with adjusted width, height, and DPI
    plt.savefig('auc_plot_288.png', dpi=300, bbox_inches='tight')

    # Display the plot
    plt.show()
    
plot_auc_curve(y_test = y_test, y_pred_proba = y_pred_proba)
