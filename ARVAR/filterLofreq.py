import pandas as pd

refDf=pd.read_csv('GA-EHC-2884X_L1_bbmap-1_pos-filter (1).txt', sep = '\t')

resDf=pd.read_csv('GA-EHC-2884X_L1_bbmap-1_lofreq-output.txt', sep = '\t')

resDf['Ref_Al_RelPos'] = 'NaN'
resDf['Var_Al_Relpos'] = 'NaN'


for i in range(len(resDf.index)):
  refAl = resDf['REF-NT'][i]
  varAl = resDf['VAR-NT'][i]
  # find real position of a variant nucleotide
  if len(refAl) > len(varAl):
    alPos = resDf['POSITION'][i] + 1
  else:
    alPos = resDf['POSITION'][i]
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


