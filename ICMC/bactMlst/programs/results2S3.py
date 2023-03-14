import subprocess
import os
import glob

fullDir = os.getcwd()

dirList = fullDir.split('/')

client_name = dirList[5] + "/"

runName = dirList[6] + "/"

cmd_str = "aws s3 cp --recursive ./custom_output/ s3://transfer-files-emory/ICMC/" + client_name + runName + "custom_output/"

subprocess.run(cmd_str, shell=True)

# transfer bactopia main output

def sampleNames():
  filesList = [os.path.basename(x) for x in glob.glob('./input/*R1*')]
  samplesList = []
  for i in filesList:
    sample = i.replace("_R1.fastq.gz", "") + "/"
    samplesList.append(sample)
  return samplesList

samplesList = sampleNames()

def transferBactOut():
  for sample in samplesList:
    cmd_str = "aws s3 cp --recursive ./bactopia_output/" + sample + " s3://transfer-files-emory/ICMC/" + client_name + runName + "original_bactopia_output/" + sample
    subprocess.run(cmd_str, shell=True)
    
transferBactOut()
