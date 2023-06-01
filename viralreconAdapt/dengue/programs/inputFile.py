import os
import subprocess
import pandas as pd
import glob

suffix = "_001"

def getSamplesName(suffix, path):
  filesList = glob.glob(f"{path}/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace(f"{path}/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

def inputDf(samplesNames, suffix, path):
  allR1 = []
  allR2 = []
  for i in samplesNames:
    r1 = f'{path}/{i}_R1{suffix}.fastq.gz'
    r2 = f'{path}/{i}_R2{suffix}.fastq.gz'
    allR1.append(r1)
    allR2.append(r2)
  inputDat = pd.DataFrame({"sample":samplesNames, "fastq_1":allR1, "fastq_2":allR2})
  return inputDat

samplesNames = getSamplesName(suffix = suffix, path = "input")
inputDat = inputDf(samplesNames = samplesNames, suffix = suffix, path = "input")


if len(pd.unique(inputDat['sample'])) == len(inputDat.index):
  inputDat.to_csv(path_or_buf = "input.csv", index = False)
  print("input written")
  
# water
samplesNames = getSamplesName(suffix = suffix, path = "water_input")
inputDat = inputDf(samplesNames = samplesNames, suffix = suffix, path = "water_input")
inputDat.to_csv(path_or_buf = "inputWater.csv", index = False)
