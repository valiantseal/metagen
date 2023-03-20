import subprocess
import os
import glob

fullDir = os.getcwd()
dirList = fullDir.split('/')
run_name = dirList[6] + "/"
proj_name = dirList[5] + "/"


outputFiles = glob.glob("*output*")
outputFiles.append("all_summaries")

def toS3():
  for i in outputFiles:
    curDir = i + "/"
    cmd_str = "aws s3 cp --recursive " + curDir + " s3://transfer-files-emory/Viralrecon/" + proj_name + run_name + curDir
    subprocess.run(cmd_str, shell=True)

toS3()
