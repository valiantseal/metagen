import os
import subprocess
import pandas as pd
import glob

suffix = "_001"
controlList = ["KJ-I2_S20", "KJ-I1_S19"]

def getSamplesName(suffix, dirName):
  filesList = glob.glob(f"{dirName}/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace(f"{dirName}/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

def inputDf(samplesNames, suffix, controlList, dirName):
  allR1 = []
  allR2 = []
  allControl = []
  allAnti = []
  for i in samplesNames:
    r1 = f'{dirName}/{i}_R1{suffix}.fastq.gz'
    r2 = ""
    allR1.append(r1)
    allR2.append(r2)
    if i in controlList:
      controlName = "WT_Input"
    else:
      controlName = ""
    allControl.append(controlName)
    antiB = ""
    allAnti.append(antiB)
  inputDat = pd.DataFrame({"sample":samplesNames, "fastq_1":allR1, "fastq_2":allR2, "antibody":allAnti, "control":allControl})
  return inputDat

samplesNames = getSamplesName(suffix = suffix, dirName = 'RAW_data')

inputDat = inputDf(samplesNames = samplesNames, suffix = suffix, controlList = controlList, dirName = 'RAW_data')


if len(pd.unique(inputDat['sample'])) == len(inputDat.index):
  inputDat.to_csv(path_or_buf = "input.csv", index = False)
