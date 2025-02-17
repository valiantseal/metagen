import subprocess
import os
import glob
import shutil
from multiprocessing import Pool
import time

#after prepprocess took 106m55

#61 metaseq samples 472min

t = 23
tSmall = 90

make_annotation_file = False

run_preprocess = True

# make reference annotation file
def runMakeAnnotDf(make_annotation_file):
  if make_annotation_file == True:
    cmd_str = "python programs/src/makeAnnotDfCLT.py -i ./programs/data/MN908947.3.gb"
    subprocess.run(cmd_str, shell = True)

runMakeAnnotDf(make_annotation_file = make_annotation_file)

# run preprocess   
def runPreprocess(run_preprocess):
  if run_preprocess == True:
    cmd_str = "python programs/src/preprocess.py"
    subprocess.run(cmd_str, shell = True)

runPreprocess(run_preprocess = run_preprocess)

## get iSNV files
filesList = glob.glob("input/*.sam")

def prepFiles():
  samplesList = []
  filesList = glob.glob("input/*.sam")
  for file in filesList:
    sampleName = file.replace("input/", "").replace(".sam", "")
    targDir = "process/" + sampleName + "/"
    samplesList.append(sampleName)
    os.makedirs(targDir, exist_ok = True)
    shutil.copy(file, targDir)
  return samplesList


def getSampNames():
  samplesList = []
  filesList = glob.glob("input/*.sam")
  for file in filesList:
    sampleName = file.replace("input/", "").replace(".sam", "")
    targDir = "process/" + sampleName + "/"
    samplesList.append(sampleName)
  return samplesList

samplesList = getSampNames()

def copyFiles(sampleName):
  curFile = "input/" + sampleName + ".sam"
  targDir = "process/" + sampleName + "/"
  os.makedirs(targDir, exist_ok = True)
  shutil.copy(curFile, targDir)
  
with Pool(tSmall) as pool:
  pool.map(copyFiles, samplesList)


def runGetVarFiles(sampleName):
  targDir = "process/" + sampleName + "/"
  os.chdir(targDir)
  sample = sampleName + ".sam"
  cmd_str = f'python ../../programs/src/getVarFilesCLT.py -i {sample} --mapqual 5 --basequal 5 --allele_frequency 0 \
-r ../../programs/data/MN908947.3.fna -b ../../programs/bin -t 4'
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

# run lofreq pipe
start_time = time.time()
with Pool(t) as pool:
  pool.map(runGetVarFiles, samplesList)

print("getVarFiles completed in --- %s minutes ---" % ((time.time() - start_time) /60) )

# filter and annotate lofreq results

def runFilterLofreq(sampleName):
  targDir = "process/" + sampleName + "/"
  os.chdir(targDir)
  cmd_str = "python ../../programs/src/filterLofreqCLT.py -i sample_lofreq-output.tsv \
  -r sample_pos-filter.tsv -a ../../"
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

with Pool(tSmall) as pool:
  pool.map(runFilterLofreq, samplesList)
  





  

