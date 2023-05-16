import boto3
import os
import pandas as pd
import subprocess

# get available metaseq samples
def listFiles(bucket_name, directory_prefix):
  s3 = boto3.client('s3')
  filesList = []
  response = s3.list_objects_v2(Bucket=bucket_name, Prefix=directory_prefix, Delimiter='/')
  if 'Contents' in response:
    for file in response['Contents']:
        filesList.append(file['Key'])
  return(filesList)

def formatSampleNames(samplesList, prefix):
  newList = []
  for i in samplesList:
    if "R1" in i:
      nameEdit = i.replace(prefix, "").rsplit('_', 1)[0]
      nameList = nameEdit.split("-")
      newName = "-".join(nameList[0:3])
      newList.append(newName)
  uniqList = list(pd.unique(newList))
  return uniqList

newDirList = listFiles(bucket_name = 'sarscov2-sgrnas', directory_prefix = 'vaxbt/New/')
newDirNames = formatSampleNames(samplesList = newDirList, prefix = "vaxbt/New/")

round1List = listFiles(bucket_name = 'sarscov2-sgrnas', directory_prefix = 'vaxbt/round1/')
round1Names = formatSampleNames(samplesList = round1List, prefix = "vaxbt/round1/")

round2List = listFiles(bucket_name = 'sarscov2-sgrnas', directory_prefix = 'vaxbt/round2/')
round2Names = formatSampleNames(samplesList = round2List, prefix = "vaxbt/round2/")

allNames = list(pd.unique(newDirNames + round1Names + round2Names))
allSamples = newDirList + round1List + round2List 

# get processed ampseq samples
def getAmpNames(path):
  newList = []
  df = pd.read_csv(path)
  samplesList = df["sample"].to_list()
  for i in samplesList:
    nameEdit = i.split("_", 1)[0]
    nameList = nameEdit.split("-")
    newName = "-".join(nameList[0:3])
    newList.append(newName)
  uniqList = list(pd.unique(newList))
  return uniqList

ampSamples = getAmpNames(path = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/input.csv")
    
def getOverlapSamples(list1, list2):
  newList = []
  for i in list1:
    for j in list2:
      if i in j:
        newList.append(j)
  return newList

def getOverlapNames(list1, list2):
  newList = []
  for i in list1:
    if i in list2:
      newList.append(i)
  return newList

def getNotOverlapNames(list1, list2):
  newList = []
  for i in list1:
    if i not in list2:
      newList.append(i)
  return newList

overlapSamples = getOverlapSamples(list1 = ampSamples, list2 = allSamples)

overlapNames = getOverlapNames(list1 = ampSamples, list2 = allNames)

notOverlapNames = getNotOverlapNames(list1 = ampSamples, list2 = allNames)

def downloadSamples(samplesList):
  os.makedirs('input/', exist_ok = True)
  for i in samplesList:
    samplePath = f's3://sarscov2-sgrnas/{i}'
    cmd_str = f'aws s3 cp {samplePath} ./input/'
    subprocess.run(cmd_str, shell = True)



downloadSamples(samplesList = overlapSamples)
