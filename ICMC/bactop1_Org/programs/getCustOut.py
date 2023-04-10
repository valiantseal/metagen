import os
import subprocess
import pandas as pd
import shutil

fullDir = os.getcwd()
dirList = fullDir.split('/')
runName = dirList[6]
os.makedirs("custom_output", exist_ok=True)

bactSamples = pd.read_table("fastqs.txt", sep = "\t")
samplesList = bactSamples["sample"].to_list()

# combine fastAni
def combAni(samplesList):
  combDat = pd.DataFrame()
  for sample in samplesList:
    inFile = f'fastani/bactopia-tools/fastani/fastani/{sample}/{sample}.tsv'
    sampleDf = pd.read_table(inFile, sep ="\t")
    combDat = pd.concat([combDat, sampleDf])
  return(combDat)

aniDf = combAni(samplesList)

outAni = f'./custom_output/{runName}_fastAni.csv'
aniDf.to_csv(path_or_buf = outAni, index = False, sep = ',')

