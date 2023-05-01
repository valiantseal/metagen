import subprocess
import os
import glob
import shutil
from multiprocessing import Pool
import time

t = 6

make_annotation_file = "no"

# make reference annotation file

def runMakeAnnotDf():
  cmd_str = "python programs/src/makeAnnotDfCLT.py -i ./programs/data/sequence.gb"
  subprocess.run(cmd_str, shell = True)

if make_annotation_file == "yes":
  runMakeAnnotDf()

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

samplesList = prepFiles()

def runGetVarFiles(sampleName):
  targDir = "process/" + sampleName + "/"
  os.chdir(targDir)
  sample = sampleName + ".sam"
  cmd_str = f'python ../../programs/src/getVarFilesCLT.py -i {sample} \
-r ../../programs/data/SARSCov2_Ref.fasta -b ../../programs/bin -t 4'
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

with Pool(t) as pool:
  pool.map(runFilterLofreq, samplesList)





  

