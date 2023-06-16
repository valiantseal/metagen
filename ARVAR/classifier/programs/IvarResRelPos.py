
import pandas as pd
import os
import sys

def getConsensus(metaSeq, ampSeq, minFreq, maxFreq):
  metaseq = pd.read_csv(metaSeq)
  ampseq = pd.read_csv(ampSeq)
  metaseqFilt = metaseq[(metaseq["ALT_FREQ"] >= minFreq) & (metaseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  ampseqFilt =  ampseq[(ampseq["ALT_FREQ"] >= minFreq) & (ampseq["ALT_FREQ"] <= maxFreq)].reset_index().drop(["index"], axis =1)
  targSnv = ampseqFilt["Samp_Pos_Ref_Alt"].to_list()
  ConsTest = []
  for i in range(len(metaseqFilt.index)):
    curSnv = metaseqFilt.loc[i, "Samp_Pos_Ref_Alt"]
    if curSnv in targSnv:
      curCons = 1
    else:
      curCons = 0
    ConsTest.append(curCons)
  metaseqFilt["ConsTest"] = ConsTest
  return(metaseqFilt)

# increasing lower frequency increase accuracy 
# colOption1, minFreq at least 0.02 and ampSeqIvar have the bset performance   
dfFilt = getConsensus(metaSeq = "test_consensus/metaseqIvar.csv", ampSeq = "test_consensus/ampseqConsIvar.csv", minFreq = 0.02, maxFreq = 1)

def checkSamples(path, dfFilt):
  samplesList = os.listdir(path)
  exactSamp = pd.unique(dfFilt["ExactSamp"])
  samplesList = os.listdir(path)
  print(len(exactSamp))
  print(len(samplesList))
  common_samples = [x for x in exactSamp if x in samplesList]
  print(len(common_samples))
  return common_samples

common_samples = checkSamples(path = '/home/ubuntu/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/process_par', dfFilt = dfFilt)

#
list(dfFilt)

def getCorPos(resDf, i, refAl, varAl):
  if "-" in varAl:
    alPos = resDf.loc[i, 'POS'] + 1
  else:
    alPos = resDf.loc[i, 'POS']
  return alPos

# handle deletion
def getRelPosDel(refAl, varAl, i, refSub, refSubDel):
  refRefAl = refAl + '-POS'
  Ref_Al_RelPos = refSubDel.loc[0, refRefAl]
  # find which Indel column has deletion amd record relative position
  refVarAl = len(varAl)
  if len(refSub.loc[0, 'Indel1']) == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel1-POS']
  elif len(refSub.loc[0, 'Indel2']) == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel2-POS']
  elif len(refSub.loc[0, 'Indel3']) == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel3-POS']
  else:
    Var_Al_Relpos = 'NaN'
  return Ref_Al_RelPos, Var_Al_Relpos

def getRelPosIns(refAl, varAl, i, refSub):
  refRefAl = refAl + '-POS' 
  Ref_Al_RelPos = refSub.loc[0, refRefAl]
  # find which Indel column has insertion and record relative position
  refVarAl = varAl
  if refSub.loc[0, 'Indel1'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel1-POS']
  elif refSub.loc[0, 'Indel2'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel2-POS']
  elif refSub.loc[0, 'Indel3'] == refVarAl:
    Var_Al_Relpos = refSub.loc[0, 'Indel3-POS']
  else:
    Var_Al_Relpos = 'NaN'
  return Ref_Al_RelPos, Var_Al_Relpos

# process substitutions
def getRelPosSub(refAl, varAl, i, refSub):
  refRefAl = refAl + '-POS'
  Ref_Al_RelPos = refSub.loc[0, refRefAl]
  refVarAl = varAl + '-POS'
  Var_Al_Relpos = refSub.loc[0, refVarAl]
  return Ref_Al_RelPos, Var_Al_Relpos


def getCorAndRelPos(refDf, resDf):
  corPos = []
  Ref_Al_RelPos = []
  Var_Al_RelPos = []
  for i in range(len(resDf.index)):
    refAl = resDf.loc[i, 'REF']
    varAl = resDf.loc[i, 'ALT']
    curCorPos = getCorPos(resDf = resDf, i = i, refAl = refAl, varAl = varAl)
    refSub = refDf.loc[refDf['POS'] == curCorPos].reset_index()
    if "-" in varAl:
      try:
        refSubDel = refDf.loc[refDf['POS'] == (curCorPos - 1)].reset_index()
        curRefPos, curVarPos = getRelPosDel(refAl = refAl, varAl= varAl, i = i, refSub = refSub, refSubDel = refSubDel)
      except:
        print(f'{i} result table index deletion failed')
        sys.exit(1)
    elif "+" in varAl:
      try:
        curRefPos, curVarPos = getRelPosIns(refAl = refAl, varAl= varAl, i = i, refSub = refSub)
      except:
        print(f'{i} result table index insertion failed')
        sys.exit(1)
    elif len(refAl) == len(varAl):
      try:
        curRefPos, curVarPos = getRelPosSub(refAl = refAl, varAl= varAl, i = i, refSub = refSub)
      except:
        print(f'{i} result table index substitution failed')
        sys.exit(1)
    corPos.append(curCorPos)
    Ref_Al_RelPos.append(curRefPos)
    Var_Al_RelPos.append(curVarPos)
  resDf['Position_corrected'] = corPos
  resDf['Ref_Al_RelPos'] = Ref_Al_RelPos
  resDf['Var_Al_RelPos'] = Var_Al_RelPos
  return resDf

#curSampFilt = getCorAndRelPos(refDf = curDf, resDf = curSampFilt)
#curDf = pd.read_csv("~/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/process_par/EHC-C19-3914N-L3/sample_pos-filter.tsv", sep = "\t", low_memory=False)
#bmread = pd.read_csv("~/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/process_par/EHC-C19-3914N-L3/readCountParsed.tsv", sep = "\t", low_memory=False)
#curSampFilt = dfFilt[dfFilt["ExactSamp"] == "EHC-C19-3914N-L3"].reset_index(drop = True)

def AnnotateSamples(dfFilt, common_samples, path):
  combDat = pd.DataFrame()
  for sample in common_samples:
    refPath = f'{path}/{sample}/sample_pos-filter.tsv'
    refDf = pd.read_csv(refPath, sep = "\t", low_memory=False)
    curSampFilt = dfFilt[dfFilt["ExactSamp"] == sample].reset_index(drop = True)
    curSampFilt = getCorAndRelPos(refDf = refDf, resDf = curSampFilt)
    combDat = pd.concat([combDat, curSampFilt], axis=0, ignore_index = True)
  return combDat

dfFilt = AnnotateSamples(dfFilt = dfFilt, path = "~/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/process_par", common_samples = common_samples)
dfFilt['Var_Al_RelPos'].isna().any()
dfFilt['Ref_Al_RelPos'].isna().any()

dfFilt.to_csv("Ludy_metaAmpIvar_overlapSnv_RelPos.csv", index = False)
