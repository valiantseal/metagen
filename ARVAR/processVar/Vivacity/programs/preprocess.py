import os
import subprocess
import glob
import shutil
from multiprocessing import Pool
from functools import partial

t = 2


os.makedirs("preprocess", exist_ok = True)
os.makedirs("input", exist_ok = True)

# index reference
def indRef():
  os.chdir("references")
  cmd_str = f"hisat2-build -p 16 SARSCov2_Ref.fasta SARSCov2_Ref"
  subprocess.run(cmd_str, shell = True)
  os.chdir("../")

indRef()

# 
def getSamplesName(suffix):
  filesList = glob.glob("./fastqs/*_R1*")
  samplesName = []
  for i in filesList:
    newName = i.replace("./fastqs/", "").replace(f'_R1{suffix}.fastq.gz', "")
    samplesName.append(newName)
  return samplesName

samplesNames = getSamplesName(suffix = '_001')

def splitFiles(samplesNames):
  for i in samplesNames:
    targDir = f'preprocess/{i}/'
    os.makedirs(targDir, exist_ok = True)
    curSamples = glob.glob(f"fastqs/{i}*")
    for curSample in curSamples:
      shutil.copy(curSample, targDir)

splitFiles(samplesNames)

# filter fastqs with fastp
def runFastp(suffix, sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  cmd_str = f'fastp -i {sampleName}_R1{suffix}.fastq.gz -I {sampleName}_R2{suffix}.fastq.gz \
  -o filtered_1.fastq -O filtered_2.fastq -A -l 25 -w 8 -h fastp.html'
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

# deduplication option does not work in python and/or regular bash
with Pool(t) as pool:
  suffix = '_001'
  fastp = partial(runFastp, suffix)
  pool.map(fastp, samplesNames)
  

def runHisat2(sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  cmd_str = f"hisat2 -p 8 -x ../../references/SARSCov2_Ref -1 filtered_1.fastq -2 filtered_2.fastq -S {sampleName}.sam"
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

with Pool(t) as pool:
  pool.map(runHisat2, samplesNames)
  

def copyFiles(sampleName):
  targDir = "preprocess/" + sampleName + "/"
  os.chdir(targDir)
  shutil.copy(f'{sampleName}.sam', "../../input/")
  os.chdir("../../")

with Pool(t) as pool:
  pool.map(copyFiles, samplesNames)

