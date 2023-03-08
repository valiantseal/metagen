import pandas as pd
import warnings
from pandas.core.common import SettingWithCopyWarning
import os

os.chdir('C:/Users/abomb/OneDrive - Emory University/Variant-calling-pipeline/Original_output_files_02032023')

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

# read data
refDf=pd.read_table('GA-EHC-2884X_L1_bbmap-1_pos-filter (1).txt', sep = '\t')

resDf=pd.read_table('GA-EHC-2884X_L1_bbmap-1_lofreq-output.txt', sep = '\t')

# get relative allele positions
resDf['Ref_Al_RelPos'] = 'NaN'
resDf['Var_Al_Relpos'] = 'NaN'

corPos =[]

for i in range(len(resDf.index)):
  refAl = resDf['REF-NT'][i]
  varAl = resDf['VAR-NT'][i]
  # find real position of a variant nucleotide
  if len(refAl) > len(varAl):
    alPos = resDf['POSITION'][i] + 1
  else:
    alPos = resDf['POSITION'][i]
  corPos.append(alPos)
  # select reference position
  refSub = refDf.loc[refDf['POS'] == alPos].reset_index()
  
  # find target relative position

  # handle deletion
  if len(refAl) > len(varAl):
    refRefAl = refAl[1] + '-POS' 
    resDf['Ref_Al_RelPos'][i] = refSub.loc[0, refRefAl]
    # find which indel column has deletion amd record relative position
    refVarAl = '-'  + refAl[1:]
    if refSub.loc[0, 'INDEL1'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL1-POS']
    elif refSub.loc[0, 'INDEL2'] == refVarAl:
       resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL2-POS']
    elif refSub.loc[0, 'INDEL3'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL3-POS']
    
  # process insertions
  elif len(refAl) < len(varAl):
    refRefAl = refAl + '-POS' 
    resDf['Ref_Al_RelPos'][i] = refSub.loc[0, refRefAl]
    # find which indel column has deletion amd record relative position
    refVarAl = '+'  + varAl[1:]
    if refSub.loc[0, 'INDEL1'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL1-POS']
    elif refSub.loc[0, 'INDEL2'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL2-POS']
    elif refSub.loc[0, 'INDEL3'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'INDEL3-POS']
      
  # process substitution
  elif len(refAl) == len(varAl):
    refRefAl = refAl + '-POS'
    resDf['Ref_Al_RelPos'][i] = refSub.loc[0, refRefAl]
    refVarAl = varAl + '-POS'
    resDf['Var_Al_Relpos'][i] = refSub.loc[0, refVarAl]
  
# a better way would be to create 2 lists with relative positions and add them as columns. May be redo in the future to avoid warnings

resDf['Position_corrected'] = corPos

def isConsensus(df, freq):
  consList = []
  for i in range(len(df.index)):
    if df.loc[i, 'ALLELE-FREQUENCY'] > freq:
      consList.append('consensus')
    else:
      consList.append('minority')
  df['Level'] = consList
  return(df)

resDf = isConsensus(df = resDf, freq = 0.5)

def adjFreq(df, relPos):
  newFreq = []
  allTests = []
  df['Ref_Al_RelPos'] = df['Ref_Al_RelPos'].astype(float)
  df['Var_Al_Relpos'] = df['Var_Al_Relpos'].astype(float)
  for i in range(len(df.index)):
    if df.loc[i, 'Level'] == 'consensus' and df.loc[i, 'Ref_Al_RelPos'] > relPos:
      freq = df.loc[i, 'ALLELE-FREQUENCY']
      posTest = 'Pass'
    elif df.loc[i, 'Level'] == 'consensus' and df.loc[i, 'Ref_Al_RelPos'] < relPos:
      freq = 1
      posTest = 'Fail'
    elif df.loc[i, 'Level'] == 'minority' and df.loc[i, 'Var_Al_Relpos'] > relPos:
      freq = df.loc[i, 'ALLELE-FREQUENCY']
      posTest = 'Pass'
    elif df.loc[i, 'Level'] == 'minority' and df.loc[i, 'Var_Al_Relpos'] < relPos:
      freq = 0
      posTest = 'Fail'
    newFreq.append(freq)
    allTests.append(posTest)
  df['Position_test'] = allTests
  df['Freq_adj'] = newFreq
  return(df)

resDf = adjFreq(df = resDf, relPos = 0.2)

def typeAllele(df):
  newType = []
  corRefAl = []
  corVarAl = []
  allNuclCh =[]
  for i in range(len(df.index)):
    refAl = df.loc[i, 'REF-NT']
    varAl = df.loc[i, 'VAR-NT']
    if len(refAl) > len(varAl):
      typeAl = 'Deletion'
      newRefAl = refAl[1] # confirm that reference is always a second character
      newVarAl = '-' + refAl[1:]
      nuclCh = str(df.loc[i, 'Position_corrected']) + newVarAl
    elif len(refAl) < len(varAl):
      typeAl = 'Insertion'
      newRefAl = refAl 
      newVarAl = '+' + varAl[1:]
      nuclCh = str(df.loc[i, 'Position_corrected']) + newVarAl
    elif len(refAl) == len(varAl):
      typeAl = 'Substitution'
      newRefAl = refAl 
      newVarAl = varAl
      nuclCh = refAl + str(df.loc[i, 'Position_corrected']) + varAl
    newType.append(typeAl)
    corRefAl.append(newRefAl)
    corVarAl.append(newVarAl)
    allNuclCh.append(nuclCh)
  df['Type'] = newType
  df['REF-allele_corrected'] = corRefAl
  df['VAR-allele_corrected'] = corVarAl
  df['Nucleotide_Change'] = allNuclCh
  return df

resDf = typeAllele(df = resDf)

# annotations 
annot = pd.read_csv('annotations.csv')
codons = pd.read_csv('codons_table.csv', names = ['full_name', 'aa', 'codon', 'aa2', 'full_name2'])

# annotate substitutions
def annotSubst(pos, codNumb, region, corVarAl):
  # subset codon table
  selAnot = annot.loc[(annot['Region/Gene'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  # find index of the nucleotide to change 
  chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
  # change nucleotide
  selAnot.loc[chInd, 'Ref'] = corVarAl
  # write new codon
  codList = selAnot['Ref'].tolist()
  newCod = ''.join(codList)
  newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
  if newAA == oldAA:
    mType = 'Synonimous'
  else:
    mType = 'Nonsynonimous'
  aaChange = oldAA + codNumb + newAA
  return [mType, aaChange]

# annotate insertions
def annotInsert(pos, codNumb, region, corVarAl):
  nucleotides = corVarAl[1:]
  # subset codon table
  selAnot = annot.loc[(annot['Region/Gene'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  # find index of the nucleotide to change 
  chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
  # get reference amino acid
  oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
  # find new codon and aminoacid
  if chInd == 2:
    newAA = oldAA
  elif chInd == 1:
    codList = selAnot.loc[0:1, 'Ref'].tolist()
    codList.append(nucleotides[0])
    newCod = ''.join(codList)
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  elif chInd == 0:
    codList = list(selAnot.loc[0, 'Ref'])
    codList.append(nucleotides[0:2])
    newCod = ''.join(codList)
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  if len(nucleotides) % 3 == 0:
    mType = 'Inframe'
    addNumb = '+' + str(int(len(nucleotides)/3)) + 'AA'
  else:
    mType = 'Frameshift'
    addNumb = '+' + str(int(len(nucleotides))) + 'nt'
  aaChange = ''.join([oldAA, codNumb, newAA, addNumb])
  return [mType, aaChange]

# annotate deletions
def annotDel(pos, codNumb, region, corVarAl):
  nucleotides = corVarAl[1:]
  

def annotate(df):
  allRegions = []
  allMutTypes = []
  allAaCh = []
  for i in range(len(df)):
    pos = df.loc[i, 'Position_corrected']
    region = annot.loc[annot['NT'] == pos, 'Region/Gene'].iloc[0]
    codNumb = annot.loc[annot['NT'] == pos, 'Codon#'].iloc[0]
    corVarAl = df.loc[i, 'VAR-allele_corrected']
    if region == "5'UTR" or region == "3'UTR":
      mutType = 'Synonymous'
      aaCh = region
    else:
      # handle substitutions
      if df.loc[i, 'Type'] == 'Substitution':
        resList = annotSubst(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
      # handle insertions  
      elif df.loc[ind, 'Type'] == 'Insertion':
        resList = annotInsert(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
      elif df.loc[ind, 'Type'] == 'Deletion':
        resList = annotDel(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
    
  allRegions.append(region)
  
    
