import os
import pandas as pd
import glob

os.chdir("C:/Users/abomb/OneDrive - Emory University/virus/iSNVs")

filesList = glob.glob("pbnas_ampseq_paths/*")

def combData(filesList):
    combDat = pd.DataFrame()
    for i in filesList:
        curDf = pd.read_csv(i)
        combDat = pd.concat([combDat, curDf], ignore_index=True)
    combDat.to_csv(path_or_buf="ampseq_combined_paths.csv", index=False)

combData(filesList)