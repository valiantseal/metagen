import pandas as pd
import os
import re

#os.chdir('/home/ubuntu/extraVol/Copyback/test_builds_2885y')

# read virema summary
def LastNlines(fname, N):
  resultsList = []
  with open(fname) as file:
    for line in (file.readlines() [-N:]):
      resultsList.append(line)
  return(resultsList)

# extract numbers from sentences
def extrNumb(x):
  numbList = []
  for i in range(len(x)):
    curNumb = re.findall(r'\d+', x[i])
    numbList.extend(curNumb)
  return numbList

colNames = ['total_reads', 'sing_map_reads', 'sing_map_reads_lts', 'seed', 'SF_recomb', 'vir_recomb', 'host_recomb', 'vh_recomb', 'microindels', 'threshold', 
'unIdent_insert', 'subsitution', 'unk_recomb', 'unmapped', 'time']

# run through chunks
def getAllChunks(sampList):
    finRes = pd.DataFrame()
    for chunk in sampList:
        inFile = './splitSeq2/' + chunk + '/final_virema-results.txt'
        resultsList = LastNlines(fname = inFile, N = 10)
        numbList = extrNumb(x = resultsList)
         # convert numbers to int
        intList = list(map(int, numbList))
        # make a data frame
        virRes = pd.DataFrame([intList], columns=colNames)
        virRes.insert(loc=0, column='chunk', value=chunk)
        # combine results
        finRes = pd.concat([finRes,virRes], ignore_index=True)
    return finRes
  

samplesList = os.listdir('gnuPar')

allSamples = pd.DataFrame()

# run for all samples
for sample in samplesList:
  targDir = 'gnuPar/' + sample
  os.chdir(targDir)
  chunksList = os.listdir('splitSeq2')
  combChunks = getAllChunks(sampList = chunksList).drop(['chunk'], axis=1)
  # sum numbers for all chuncks, seed, treshhold, and time would be wrong
  sumPerSamp = combChunks.sum().to_frame().transpose()
  sumPerSamp.insert(loc=0, column='sample', value=sample)
  allSamples = pd.concat([allSamples, sumPerSamp], ignore_index = True)
  os.chdir('../../')


allSamples.to_csv(path_or_buf = 'virema2SplitRes.csv', index = False)
