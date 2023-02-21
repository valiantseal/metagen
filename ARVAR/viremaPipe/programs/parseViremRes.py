import pandas as pd
import os
import re

os.chdir('/home/ubuntu/extraVol/Copyback/test_builds_2885y')

samplesList = os.listdir('gnuPar')

# read
def LastNlines(fname, N):
  resultsList = []
  with open(fname) as file:
    for line in (file.readlines() [-N:]):
      resultsList.append(line)
  return(resultsList)

def extrNumb(x):
  numbList = []
  for i in range(len(x)):
    curNumb = re.findall(r'\d+', x[i])
    numbList.extend(curNumb)
  return numbList

colNames = ['total_reads', 'sing_map_reads', 'sing_map_reads_lts', 'seed', 'SF_recomb', 'vir_recomb', 'host_recomb', 'vh_recomb', 'microindels', 'threshold', 
'unIdent_insert', 'subsitution', 'unk_recomb', 'unmapped', 'time']

# run through samples
def getAllSamples(sampList):
    finRes = pd.DataFrame()
    for sample in sampList:
        inFile = './gnuPar/' + sample + '/final_virema-results.txt'
        resultsList = LastNlines(fname = inFile, N = 10)
        numbList = extrNumb(x = resultsList)
        intList = list(map(int, numbList))
        virRes = pd.DataFrame([intList], columns=colNames)
        virRes.insert(loc=0, column='sample', value=sample)
        finRes = pd.concat([finRes,virRes], ignore_index=True)
    return finRes
  
finalRes = getAllSamples(sampList = samplesList)

