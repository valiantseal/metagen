import subprocess
import os
import glob
import shutil
from multiprocessing import Pool
import time

start_time = time.time()

#after prepprocess took 106m55

t = 12
tSmall = 46

make_annotation_file = "true"

run_preprocess = "true"

install_packages = "true"

# run package instalations
def installPackages(install_packages):
  if install_packages == "true":
    cmd_str = "python programs/src/installPackages.py"
    subprocess.run(cmd_str, shell = True)

installPackages(install_packages = install_packages)

# make reference annotation file
def runMakeAnnotDf(make_annotation_file):
  if make_annotation_file == "true":
    cmd_str = "python programs/src/makeAnnotDfCLT.py -i ./programs/data/sequence.gb"
    subprocess.run(cmd_str, shell = True)

runMakeAnnotDf(make_annotation_file = make_annotation_file)

# run preprocess   
def runPreprocess(run_preprocess):
  if run_preprocess == "true":
    cmd_str = "python programs/src/preprocess.py"
    subprocess.run(cmd_str, shell = True)

runPreprocess(run_preprocess = run_preprocess)

print("preprocessing completed in --- %s minutes ---" % ((time.time() - start_time) /60) )

# get iSNV files
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
  shutil.move(curFile, targDir)
  
with Pool(tSmall) as pool:
  pool.map(copyFiles, samplesList)


def runGetVarFiles(sampleName):
  targDir = "process/" + sampleName + "/"
  os.chdir(targDir)
  sample = sampleName + ".sam"
  cmd_str = f'python ../../programs/src/getVarFilesCLT.py -i {sample} --mapqual 5 --basequal 5 --allele_frequency 0 \
-r ../../programs/data/SARSCov2_Ref.fasta -b ../../programs/bin -t 4'
  subprocess.run(cmd_str, shell = True)
  os.chdir("../../")

# run lofreq pipe
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
  

print("Whole pipeline completed in --- %s minutes ---" % ((time.time() - start_time) /60) )



  

