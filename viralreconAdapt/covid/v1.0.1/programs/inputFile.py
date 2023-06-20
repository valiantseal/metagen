import os
import subprocess
import pandas as pd
import glob

suffix = "_001"

def getSamplesName(suffix):
  filesList = glob.glob("./input/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace("./input/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

def inputDf(samplesNames, suffix):
  allR1 = []
  allR2 = []
  for i in samplesNames:
    r1 = f'input/{i}_R1{suffix}.fastq.gz'
    r2 = f'input/{i}_R2{suffix}.fastq.gz'
    allR1.append(r1)
    allR2.append(r2)
  inputDat = pd.DataFrame({"sample":samplesNames, "fastq_1":allR1, "fastq_2":allR2})
  return inputDat

samplesNames = getSamplesName(suffix = suffix)

inputDat = inputDf(samplesNames = samplesNames, suffix = suffix)


if len(pd.unique(inputDat['sample'])) == len(inputDat.index):
  inputDat.to_csv(path_or_buf = "input.csv", index = False)



