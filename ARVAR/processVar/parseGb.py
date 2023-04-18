import time
start_time = time.time()

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

def findUtrInd(genFile):
  for i in range(len(genFile)):
    if "5'UTR" in genFile[i]:
      seqStart = i
    if "3'UTR" in genFile[i]:
      seqEnd = i +1
  return seqStart, seqEnd

seqStart, seqEnd = findUtrInd(genFile)


def ParseFile(seqStart, seqEnd):
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

parsedFile = ParseFile(seqStart, seqEnd)

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
#max(266,806) < min(21555,2719)

def findSeq(genFile):
  for i in range(len(genFile)):
    if "ORIGIN" in genFile[i]:
      seqStart = i + 1
    if "//" in genFile[i]:
      seqEnd = i
  return seqStart, seqEnd

fastStart, fastEnd = findSeq(genFile)

def getFasta(fastStart, fastEnd):
  results = ""
  for i in range(fastStart, fastEnd):
    newLine = genFile[i].rstrip()
    editLine = ' '.join( newLine.split()).upper()
    lineList =  editLine.split(" ")
    filtLine = ''.join(lineList[1:])
    results = results + filtLine
    
  return results

fastaSeq = getFasta(fastStart, fastEnd)

fastaDf = pd.DataFrame(list(fastaSeq), columns =['Ref'])
fastaDf.insert(loc = 0, column = "NT", value = fastaDf.index +1)

def splitFeatures(parsedResult):
  combDf = pd.DataFrame()
  for i in range(len(parsedResult.index)):
    newDf = pd.DataFrame(index=pd.Series(range(parsedResult.loc[i, "Start_pos"],parsedResult.loc[i, "End_pos"] +1), dtype = 'float64'))
    newDf["NT"] = newDf.index
    newDf["Title"] = parsedResult.loc[i, "Title"]
    newDf["Gene"] = parsedResult.loc[i, "Gene"]
    newDf["Product"] = parsedResult.loc[i, "Product"]
    newDf["Protein_id"] = parsedResult.loc[i, "Protein_id"]
    newDf["Length"] = parsedResult.loc[i, "Length"]
    combDf = pd.concat([combDf, newDf], ignore_index=True)
  sortDf = combDf.sort_values(by=['NT']).reset_index().drop(['index'], axis = 1)
  return(sortDf)

refDf = splitFeatures(parsedResult)

def filterNA(x):
  newList = []
  for i in x:
    if i != 'NaN':
      newList.append(i)
  return newList

def getMinVar(positionList, header, nuclRef):
  for position in positionList:
      minRef = nuclRef[nuclRef["Length"] == position]
      combVar = list(np.unique(minRef[header]))
      if any(val != "NaN" for val in combVar):
        combVar = filterNA(combVar)
        combStr = ",".join(combVar)
        break
  if 'combStr' not in locals():
    combStr = "NaN"
  return combStr

def annotateNucl(refDf, i):
  nuclRef = refDf[refDf["NT"] == i]
  if len(nuclRef.index) > 0:
    nuclRef = nuclRef.sort_values(by = ["Length"])
    #minRef = nuclRef[nuclRef["Length"] == min(nuclRef['Length'])] 
    headList = ["Title", "Gene", "Product", "Protein_id"]
    positionList = list(np.unique(nuclRef["Length"]))
    refList = []
    minList = []
    for header in headList:
      combVar = np.unique(nuclRef[header])
      combStr = ",".join(filterNA(combVar))
      if combStr == '':
        combStr = "NaN"
      refList.append(combStr)
    for header in headList:
      minVar = getMinVar(positionList, header, nuclRef)
      minList.append(minVar)
  else:
    refList = ["NCR"] *4
    minList = ["NCR"] *4
  return refList, minList

def annotFastaDf(refDf, fastaDf):
  nuclPositions = fastaDf["NT"].to_list()
  allFullRef = []
  allMinRef = []
  for i in nuclPositions:
    fullRef, minRef = annotateNucl(refDf, i)
    allFullRef.append(fullRef)
    allMinRef.append(minRef)
  
  fullRefDf = pd.DataFrame(allFullRef, columns = ["All_Titles", "All_Genes", "All_Products", "All_Protein_IDs"])
  minRefDf = pd.DataFrame(allMinRef, columns = ["Min_Titles", "Min_Genes", "Min_Products", "Min_Protein_IDs"])
  
  annotDf = pd.concat([fastaDf, minRefDf, fullRefDf], axis = "columns")
  return(annotDf)

annotDf = annotFastaDf(refDf, fastaDf)

def prevElem(my_list, i, n):
  for elem in range(i-1, -1, -1):
    prevElem = my_list[elem]
    if prevElem != n:
      newElem = prevElem
      break
  return newElem

def nextElem(my_list, i, n):
  for elem in range(i, len(my_list)):
    nextElem = my_list[elem]
    if nextElem != n:
      newElem = nextElem
      break
  return newElem

def reformPrevNext(my_list, n):
  resList = []
  for i in range(len(my_list)):
    curVal = my_list[i]
    if curVal == n:
      prevVal = prevElem(my_list, i, n)
      nextVal = nextElem(my_list, i, n)
      newVal = "_".join([str(prevVal), str(curVal), str(nextVal)])
    else:
      newVal = str(curVal)
    resList.append(newVal)
  return resList

annotDf["All_Products_Edit"] = reformPrevNext(my_list = annotDf["All_Products"].to_list(), n = "NCR")

#q1 = annotDf[annotDf["Min_Products"] == "RNA-dependent RNA polymerase"]
  
def findStartCod(annotDf):
  for i in range(len(annotDf.index)):
    if "5'UTR" not in annotDf.loc[i, "All_Titles"]:
      targNuc = annotDf.loc[i, "NT"]
      break
  return targNuc

def findEndCod(annotDf):
  for i in range(len(annotDf.index)):
    if "3'UTR" in annotDf.loc[i, "All_Titles"]:
      targNuc = annotDf.loc[i-1, "NT"]
      break
  return targNuc

startCod = findStartCod(annotDf)
    
endCod = findEndCod(annotDf)

def addCodons(fastaSeq, startCod, endCod, annotDf ):
  utr5Seq = ["NCR"] * (startCod-1)
  utr3Seq = (len(fastaSeq) - endCod) * ["NCR"]
  subStr = fastaSeq[(startCod -1):endCod]
  # split string in a list of 3 elements
  substrList = [subStr[i:i+3] for i in range(0, len(subStr), 3)]
  cdsList = []
  for item in substrList:
    cdsList.extend([item]*3)
  seqList = utr5Seq + cdsList + utr3Seq
  len(seqList) == len(fastaSeq)
  annotDf["Ref_Codon"] = seqList
  # codons for coding sequence
  #posList = list(range(1, 1 + int((len(subStr)/3))))
  #allPos = []
  #for item in posList:
    #allPos.extend([item]*3)
  #codNumb = utr5Seq + allPos + utr3Seq
  #annotDf["Codon#"] = codNumb
  return(annotDf)

annotDf = addCodons(fastaSeq, startCod, endCod, annotDf)

#annotDf['Region/Gene'] = annotDf['Min_Genes'] + '_' + annotDf['Min_Products']

def addCodNumbPerBase(annotDf):
  editDf = pd.DataFrame()
  productsList = list(pd.unique(annotDf["All_Products_Edit"]))
  for i in productsList:
    curProd = annotDf[annotDf['All_Products'] == i].reset_index().drop(['index'], axis =1)
    curProd["Codon#"] = "NaN"
    for j in range(len(curProd.index)):
      if (curProd.loc[j, "All_Products"] == "NaN") or (curProd.loc[j, "All_Products"] == "NCR"):
        curProd.loc[j, "Codon#"] = "NCR"
      else:
        if j < 3:
          curProd.loc[j, "Codon#"] = 1
        else:
          newCod = curProd.loc[j-3, "Codon#"] + 1
          curProd.loc[j, "Codon#"] = newCod
        
    editDf = pd.concat([editDf, curProd])
  editDf = editDf.sort_values(by = ["NT"]).reset_index().drop(["index"], axis =1)
  return(editDf)

annotDf = addCodNumbPerBase(annotDf)

def transLateCod(annotDf):
  aaList = []
  codons = pd.read_csv('codons_table.csv', names = ['full_name', 'aa', 'codon', 'aa2', 'full_name2'])
  for i in range(len(annotDf.index)):
    codon = annotDf.loc[i, "Ref_Codon"]
    if codon == "NCR" or codon == "NaN":
      newAA = "NCR"
    else:
      newAA = codons.loc[codons['codon'] == codon, 'aa'].iloc[0]
    aaList.append(newAA)
  annotDf["Ref_AA"] = aaList
  return(annotDf)

annotDf = transLateCod(annotDf)

annotDf = annotDf.rename(columns = {'Min_Products' : 'Region/Gene'})
annotDf.to_csv(path_or_buf = "py_annotations.csv", index = False)
    
print("--- %s minutes ---" % ((time.time() - start_time) /60) )
