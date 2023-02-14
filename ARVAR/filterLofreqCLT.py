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
from pandas.core.common import SettingWithCopyWarning

"""
a better way would be to create 2 lists with relative positions and add them as columns. 
May be redo in the future to avoid warnings

"""
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
parser=argparse.ArgumentParser(description='format metadata for Nextrain')

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
  "--filter",
  type=float,
  nargs="?",
  default=None,
  help="filter resullts to include relative variant position higher or equal to INT"
)

# parse arguments
args = parser.parse_args()

df = args.input_path
ref = args.reference
output = args.output_path
filt = args.filter

# main function
refDf=pd.read_table(ref, sep = '\t')

resDf=pd.read_table(df, sep = '\t')

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

# filter and write
if filt == None:
  finalDf = resDf
else:
  filtNa = resDf.loc[resDf['Var_Al_Relpos'] != 'NaN']
  filtNa['Var_Al_Relpos'] = filtNa['Var_Al_Relpos'].astype(float)
  finalDf = filtNa.loc[filtNa['Var_Al_Relpos'] >= filt]

finalDf.to_csv(path_or_buf= output, index=False)
