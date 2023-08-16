import shutil
import pandas as pd
import os
from multiprocessing import Pool, freeze_support

os.chdir("C:/Users/abomb/OneDrive - Emory University/virus/iSNVs")

def copyFiles(df, outDir):
    os.makedirs(outDir, exist_ok=True)
    curDf = pd.read_csv(df)
    for i in range(len(curDf.index)):
        inFile = curDf.loc[i, 'Pbnas_path']
        outFile = outDir + "/" + os.path.basename(inFile)
        shutil.copy(inFile, outFile)

def getPaths(df, outDir):
    os.makedirs(outDir, exist_ok=True)
    combPaths = []
    curDf = pd.read_csv(df)
    for i in range(len(curDf.index)):
        inFile = curDf.loc[i, 'Pbnas_path']
        combPaths.append(inFile)
    return combPaths

def copyFastqPar(inFile):
    outFile = "amp_fail_files/" + os.path.basename(inFile)
    shutil.copy(inFile, outFile)

curPaths = getPaths(df = "amp_uniqueFailFiltr.csv", outDir = "amp_fail_files")

if __name__ == "__main__":
    #copyFiles(df = "amp_uniqueMiss.csv", outDir = "amp_miss_files")
    freeze_support()
    with Pool(12) as pool:
        pool.map(copyFastqPar, curPaths)