# check if the required packages can be found
import sys
import subprocess
import pkg_resources

required = {'pandas', 'argparse'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

# install required packages
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)

### main program
import pandas as pd
import argparse
import warnings
import math

"""
a better way would be to create 2 lists with relative positions and add them as columns. 
May be redo in the future to avoid warnings

"""

parser=argparse.ArgumentParser(description='format and filter Lofreq output')

# add arguments
parser.add_argument(
  "-i",
  "--input_path",
  type=str,
  help="Path to the input file, required"
)

parser.add_argument(
  "-r",
  "--reference",
  type=str,
  help="path to a reference table"
)

parser.add_argument(
  "-o",
  "--output_path",
  type=str,
  nargs="?",
  default="filtered.csv",
  help="path to output file"
)

parser.add_argument(
  "-f",
  "--filter_res",
  type=float,
  nargs="?",
  default=None,
  help="filter resullts to include only relative variant position higher or equal to Float"
)

parser.add_argument(
  "-c",
  "--consensus",
  type=float,
  nargs="?",
  default=0.5,
  help="variant frequencies higher than a value will be considered consensus variants, Float"
)

parser.add_argument(
  "-p",
  "--relative_pos",
  type=float,
  nargs="?",
  default=0.2,
  help="relative position threshold for adjusted frequency, Float"
)

parser.add_argument(
  "-s",
  "--split",
  type=str,
  nargs="?",
  default='F',
  help="write variants that failed the position test separately, string T or F"
)

parser.add_argument(
  "-a",
  "--annot_dir",
  type=str,
  nargs="?",
  default='.',
  help="path to the directory with annotation files, default '.' "
)

# parse arguments
args = parser.parse_args()

df = args.input_path
ref = args.reference
output = args.output_path
filt = args.filter_res
cons = args.consensus
rel = args.relative_pos
splitVar = args.split
annot_dir = args.annot_dir
## main function

# read data
refDf=pd.read_table(ref, sep = '\t')

resDf=pd.read_table(df, sep = '\t')

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
    # find which Indel column has deletion amd record relative position
    refVarAl = '-'  + refAl[1:]
    if refSub.loc[0, 'Indel1'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel1-POS']
    elif refSub.loc[0, 'Indel2'] == refVarAl:
       resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel2-POS']
    elif refSub.loc[0, 'Indel3'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel3-POS']
    
  # process insertions
  elif len(refAl) < len(varAl):
    refRefAl = refAl + '-POS' 
    resDf['Ref_Al_RelPos'][i] = refSub.loc[0, refRefAl]
    # find which Indel column has deletion amd record relative position
    refVarAl = '+'  + varAl[1:]
    if refSub.loc[0, 'Indel1'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel1-POS']
    elif refSub.loc[0, 'Indel2'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel2-POS']
    elif refSub.loc[0, 'Indel3'] == refVarAl:
      resDf['Var_Al_Relpos'][i] = refSub.loc[0, 'Indel3-POS']
      
  # process substitution
  elif len(refAl) == len(varAl):
    refRefAl = refAl + '-POS'
    resDf['Ref_Al_RelPos'][i] = refSub.loc[0, refRefAl]
    refVarAl = varAl + '-POS'
    resDf['Var_Al_Relpos'][i] = refSub.loc[0, refVarAl]
  
# a better way would be to create 2 lists with relative positions and add them as columns. May be redo in the future to avoid warnings
resDf['Position_corrected'] = corPos

# find if the variant is a consensus level
def isConsensus(df, freq):
  consList = []
  for i in range(len(df.index)):
    if df.loc[i, 'ALLELE-FREQUENCY'] > freq:
      consList.append('consensus')
    else:
      consList.append('minority')
  df['Level'] = consList
  return(df)

resDf = isConsensus(df = resDf, freq = cons)

# adjust variants frequency based on threshold
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

# correct reference and variant alleles, add nucleotide change
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

## annotations 
# do we need these files as user input options?
annot = pd.read_csv(f'{annot_dir}/py_annotations.csv')
codons = pd.read_csv(f'{annot_dir}/programs/data/codons_table.csv', names = ['full_name', 'aa', 'codon', 'aa2', 'full_name2'])

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
    mType = 'Synonymous'
  else:
    mType = 'Nonsynonymous'
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
    codList = selAnot['Ref'].tolist()
    codList.insert(2, nucleotides)
    newSeq = ''.join(codList)
    newCod = newSeq[0:3]
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  elif chInd == 0:
    codList = selAnot['Ref'].tolist()
    codList.insert(1, nucleotides)
    newSeq = ''.join(codList)
    newCod = newSeq[0:3]
    newAA = codons.loc[codons['codon'] == newCod, 'aa'].iloc[0]
  # inframe or not  
  if len(nucleotides) % 3 == 0:
    mType = 'Inframe'
    addNumb = '+' + str(int(len(nucleotides)/3)) + 'AA'
  else:
    mType = 'Frameshift'
    addNumb = '+' + str(int(len(nucleotides))) + 'nt'
  # if to write a new aa  
  if chInd == 2:
    aaChange = ''.join([oldAA, codNumb, addNumb])
  else:
    aaChange = ''.join([oldAA, codNumb, newAA, addNumb])
  return [mType, aaChange]


# annotate deletions
def annotDel(pos, codNumb, region, corVarAl):
  nucleotides = corVarAl[1:]
  # subset codon table
  selAnot = annot.loc[(annot['Region/Gene'] == region) & (annot['Codon#'] == str(codNumb))].reset_index().drop(['index'], axis=1)
  # find index of the nucleotide to change 
  chInd = selAnot.index[selAnot['NT'] == pos].astype('int')
  # deletions never need a new amino acid or a check at what codon position nucleotide is deleted
  # needs old amino acid, codon number, -numb of aa or nt
  oldAA = selAnot.loc[chInd, 'Ref_AA'].iloc[0]
  if len(nucleotides) % 3 == 0:
    mType = 'Inframe'
    addNumb = '-' + str(int(len(nucleotides)/3)) + 'AA'
  else:
    mType = 'Frameshift'
    addNumb = '-' + str(int(len(nucleotides))) + 'nt'
  aaChange = ''.join([oldAA, codNumb, addNumb])
  return [mType, aaChange]

# use functions to annotate the table
def annotate(df):
  allRegions = []
  allMutTypes = []
  allAaCh = []
  for i in range(len(df)):
    pos = df.loc[i, 'Position_corrected']
    region = annot.loc[annot['NT'] == pos, 'Region/Gene'].iloc[0]
    codNumb = annot.loc[annot['NT'] == pos, 'Codon#'].iloc[0]
    corVarAl = df.loc[i, 'VAR-allele_corrected']
    if region == "5'UTR" or region == "3'UTR" or region == 'E/NCR' or region == 'NaN' or pd.isna(region) or region == 'NCR':
      mutType = 'Synonymous'
      aaCh = region
    else:
      # handle substitutions
      if df.loc[i, 'Type'] == 'Substitution':
        resList = annotSubst(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
      # handle insertions  
      elif df.loc[i, 'Type'] == 'Insertion':
        resList = annotInsert(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
      # handle deletions
      elif df.loc[i, 'Type'] == 'Deletion':
        resList = annotDel(pos=pos, codNumb=codNumb, region=region, corVarAl=corVarAl)
      mutType = resList[0]
      aaCh = resList[1]
    # combine lists
    allRegions.append(region)
    allMutTypes.append(mutType)
    allAaCh.append(aaCh)
    # make new columns
  df['Region'] = allRegions
  df['Mutation_type'] = allMutTypes
  df['AA_change'] = allAaCh
  return df

resDf = annotate(df=resDf)

# add  Pi*Ln(Pi) for shannon index 
def ShPi(df):
  pi = []
  for i in range(len(df)):
    if df.loc[i, 'Freq_adj'] == 0:
      newPi = 0
    else:
      newPi = df.loc[i, 'Freq_adj'] * math.log(df.loc[i, 'Freq_adj'])
    pi.append(newPi)
  df['Pi*Ln(Pi)'] = pi
  return df

resDf = ShPi(df = resDf)

## filter and write
if filt == None:
  finalDf = resDf
else:
  filtNa = resDf.loc[resDf['Var_Al_Relpos'] != 'NaN']
  filtNa['Var_Al_Relpos'] = filtNa['Var_Al_Relpos'].astype(float)
  finalDf = filtNa.loc[filtNa['Var_Al_Relpos'] >= filt]

  
if splitVar == 'T':
  posRes = finalDf.loc[finalDf['Position_test'] == 'Pass']
  negRes = finalDf.loc[finalDf['Position_test'] == 'Fail' ]
  negRes.to_csv(path_or_buf= 'failed_position_test.csv', index=False)
  posRes.to_csv(path_or_buf= output, index=False)
else:
  finalDf.to_csv(path_or_buf= output, index=False)
