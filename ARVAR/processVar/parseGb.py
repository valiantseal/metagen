import os
import pandas as pd
import numpy as np

os.chdir('C:/Users/abomb/OneDrive - Emory University/Variant-calling-pipeline/Original_output_files_02032023')

###

def getReadCountRes(filePath):
  with open(filePath, "r") as f:
    lines=f.readlines()
  return(lines)

genFile = getReadCountRes(filePath = "sequence.gb")

def findUtrInd():
  for i in range(len(genFile)):
    if "5'UTR" in genFile[i]:
      seqStart = i
    if "3'UTR" in genFile[i]:
      seqEnd = i +1
  return seqStart, seqEnd

seqStart, seqEnd = findUtrInd()


def ParseFile():
  results = []
  for i in range(seqStart, seqEnd):
    newLine = genFile[i].split(" ")
    if newLine[5] == '':
      strLine = ' '.join(newLine)
      editLine = ' '.join(strLine.split())
    else:
      strLine = ' '.join(newLine)
      editLine = ' '.join(strLine.split())
      editLine = "__" + editLine
    results.append(editLine)

  resStr = ' '.join(results)
  editList = resStr.split("__")
  editList.pop(0)
  return(editList)

parsedFile = ParseFile()

def getHeader(curList):
    headerList = curList[0].split(" ")
    curTitle = headerList[0]
    if "join(" not in headerList[1]:
      startPos, endPos = headerList[1].split("..")
    else:
      posList = headerList[1].replace("join(", "").replace(")", "").replace("..", ",").split(",")
      startPos = posList[0]
      endPos = posList[(len(posList) - 1)]
    return [curTitle, startPos, endPos]
    

def extractInfo(parsedList):
  allResults = []
  for i in range(len(parsedList)):
    curList = parsedList[i].split(" /")
    resList = getHeader(curList = curList)
    # extract elements from the list
    for elem in curList:
      if "gene=" in elem:
        curGene = eval(elem.replace("gene=", ""))
      elif "product=" in elem:
        curProduct = eval(elem.replace("product=", ""))
      elif "protein_id=" in elem:
        curProtId = eval(elem.replace("protein_id=", ""))
      elif "note=" in elem:
        curNote = eval(elem.replace("note=", ""))
      # deal with missing data  
      if 'curGene' not in locals():
        curGene = "NaN"
      if 'curProduct' not in locals():
        curProduct = "NaN"
      if 'curProtId' not in locals():
        curProtId = "NaN"
      if "curNote" not in locals():
        curNote = "NaN"
    
    resList.extend([curGene, curProduct, curProtId, curNote])
    del(curGene, curProduct, curProtId, curNote)
    allResults.append(resList)
    resultDf = pd.DataFrame( allResults, columns = ["Title", "Start_pos", "End_pos", "Gene", "Product", "Protein_id", "Note"])
  return resultDf

parsedResult = extractInfo(parsedList = parsedFile)

#parsedResDed= parsedResult.drop_duplicates()

parsedResult["Range"] = parsedResult["Start_pos"] + "-" + parsedResult["End_pos"]
parsedResult[["Start_pos", "End_pos"]] = parsedResult[["Start_pos", "End_pos"]].apply(pd.to_numeric)  
parsedResult["Length"] = parsedResult["End_pos"] - parsedResult["Start_pos"] +1
    
###

df = pd.DataFrame(index = list(np.arange(10,20)))
df = pd.DataFrame(index=pd.Series(range(10,20)))
df["Gene"] = "Gene"
df["Position"] = df.index

df1 =  pd.DataFrame(index=pd.Series(range(15,25)))
df1["Gene"] = "Prot"
df1["Position"] = df1.index

df3 = pd.merge(df, df1, how = "outer", on = "Position")


###
max(266,806) < min(21555,2719)

# solution for each row i and j check if ranges overlap, if yes, and one range is higher it will be grand titel if not paste together, also start and end should not be identical

# make main table with all cds positions, make 3 tables cds, genes, other, if position in other and fetures not NA fill if NA fill with next one
# to avoid duplicates make column with rangees, subset by unique range if fetures do not match paste them.

smallSegments = parsedResult[(parsedResult['Title'] != 'CDS') & (parsedResult['Title'] != 'gene')]

rangesList = list(np.unique(smallSegments["Range"]))

