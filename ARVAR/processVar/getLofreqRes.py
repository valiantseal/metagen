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

def parseReadCount():
  
  for i in range(len(readCountRes)):
    
    

