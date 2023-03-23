import os
import subprocess
from itertools import islice
import pandas as pd
import numpy as np


threads = '8'
mapqual = '20'
basequal = '20'
refFasta = './references/MN908947.3.fna'

# get reference name and length from alignment file bowtie2 can print lines nicely       
def readFileLines(f):
  linesList =[]
  N = 2 
  with open(f, "r") as file:
    for i in range(N):
      line = next(file).replace('\t', ' ')
      linesList.append(line)
  targLineList = linesList[1].split(' ')
  ref = targLineList[1].replace("SN:", "")
  alLen = targLineList[2].replace("LN:", "").replace("\n", "")
  return ref, alLen

# easier to understand version 
def getRefLen(f):
  with open(f) as input_file:
    head = list(islice(input_file, 2))
  targLine = head[1].split('\t')
  ref = targLine[1].replace("SN:", "")
  alLen = targLine[2].replace("LN:", "").replace("\n", "")
  return ref, alLen

refName, refLen = getRefLen(f = "input/output.sam")

#Convert the sam file into a bam file containing only mapped reads
def sam2Bam(sam):
  cmd_str = 'samtools view -@ ' + threads + ' -bu -F 4 ' + sam + ' -o - | samtools sort -@ ' + threads + ' - -o output.bam'
  subprocess.run(["bash", "-c", cmd_str],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE,
                          text=True)
                          
sam2Bam(sam = "input/output.sam")

# index bam file
def indexBam(bam):
  cmd_str = 'samtools index -@ ' + threads + ' ' + bam
  subprocess.run(cmd_str, shell=True, stderr=subprocess.PIPE)

indexBam(bam = 'output.bam')

# bam readcount
def bamReadCount(bam):
  cmd_str = 'bam-readcount/bam-readcount -w1 -q %s -b %s -f %s %s %s:1-%s > readcount.tsv' % (mapqual, basequal, refFasta, bam, refName, refLen)
  subprocess.run(cmd_str, shell=True)

bamReadCount(bam = 'output.bam')

# parse readcount file 
def getReadCountRes():
  with open("readcount.tsv", "r") as f:
    lines=f.readlines()
    result=[]
    for x in lines:
      result.append(x.split('\t'))
  return(result)

readCountRes = getReadCountRes()

# parse results
def getColNumbers():
  colNumb =[]
  for i in readCountRes:
    colNumb.append(len(i))
  uniqNumb = list(np.unique(colNumb))
  return uniqNumb

colNumbers = getColNumbers()

# extract info for nucleotides
def parseNucleotide(i, col):
  nuclRes = readCountRes[i][col].split(':')
  selRes = list(nuclRes[1].split(' ')) + nuclRes[5:9]
  return selRes

# extract infor for indels
def parseIndels(i, col):
  emptList = ["NaN"] * 6
  if len(readCountRes[i]) > col:
    nuclRes = readCountRes[i][col].split(':')
    selRes = nuclRes[0:2] + nuclRes[5:9]
  else:
    selRes = emptList
  return selRes

# parse whole bam readcount
def parseReadCount():
  allMainRes = []
  allAres = []
  allCres = []
  allGres = []
  allTres = []
  allIndel1 = []
  for i in range(len(readCountRes)):
    mainRes = readCountRes[i][0:4]
    aRes = parseNucleotide(i = i, col = 5)
    cRes = parseNucleotide(i = i, col = 6)
    gRes = parseNucleotide(i = i, col = 7)
    tRes = parseNucleotide(i = i, col = 8)
    indel1 = parseIndels(i = i , col = 10)
    #combine lists
    allMainRes.append(mainRes)
    allAres.append(aRes)
    allCres.append(cRes)
    allGres.append(gRes)
    allTres.append(tRes)
    allIndel1.append(indel1)
    
  # make dataframes  
  mainData = pd.DataFrame(allMainRes, columns = ['chr',	'position',	'reference_base',	'depth'])
  aDf = pd.DataFrame(allAres, columns = ['A_count', 'A_num_plus_strand', 'A_num_minus_strand', 'A_avg_pos_as_fraction', 'A_avg_num_mismatches_as_fraction'])
  cDf = pd.DataFrame(allCres, columns = ['C_count', 'C_num_plus_strand', 'C_num_minus_strand', 'C_avg_pos_as_fraction', 'C_avg_num_mismatches_as_fraction'])
  gDf = pd.DataFrame(allGres, columns = ['G_count', 'G_num_plus_strand', 'G_num_minus_strand', 'G_avg_pos_as_fraction', 'G_avg_num_mismatches_as_fraction'])
  tDf = pd.DataFrame(allTres, columns = ['T_count', 'T_num_plus_strand', 'T_num_minus_strand', 'T_avg_pos_as_fraction', 'T_avg_num_mismatches_as_fraction'])
  indel1Df = pd.DataFrame(allIndel1, columns = ['Indel1_base', 'Indel1_count', 'Indel1_num_plus_strand', 'Indel1_num_minus_strand', 'Indel1_avg_pos_as_fraction', 'Indel1_avg_num_mismatches_as_fraction'])
  
  # combine tables  
  combData = pd.concat([mainData, aDf, cDf, gDf, tDf, indel1Df], axis="columns")
  return combData

mainData = parseReadCount()

#subData = mainData[['A_count', 'T_count', 'A_avg_pos_as_fraction']]
#mainData.to_csv(path_or_buf = 'readcountFormat.csv', index = False)
    
    

