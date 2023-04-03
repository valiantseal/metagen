import pandas as pd
import pytest
import os
#from itertools import islice
#import subprocess

#lofreq 2.1.5

def test_ReadCount(df1, df2):
  pyReadC = pd.read_table(df1, sep  = "\t")
  shReadC = pd.read_table(df2, sep = "\t")
  pyReadSort = pyReadC.sort_values(by = ['POS'])
  shReadSort = shReadC.sort_values(by = ['POS'])
  #assert pyReadSort.equals(shReadSort)
  if pyReadSort.equals(shReadSort):
    print("bam-readcount pass")
  else:
    print("bam-readcount FAIL")

test_ReadCount(df1 = "sample_pos-filter.tsv" , df2 = "test_rose/position-filters/GA-EHC-2884X_L1_pos-filter.txt")

###
def test_lofreqFinal(df1, df2):
  pyReadC = pd.read_table(df1, sep  = "\t")
  shReadC = pd.read_table(df2, sep = "\t")
  pyReadSort = pyReadC.sort_values(by = ['POSITION'])
  shReadSort = shReadC.sort_values(by = ['POSITION'])
  pyReadSort['Pos_Ref_Alt'] = pyReadSort['POSITION'].astype(str) + pyReadSort['REF-NT'] + pyReadSort['VAR-NT']
  shReadSort['Pos_Ref_Alt'] = shReadSort['POSITION'].astype(str) + shReadSort['REF-NT'] + shReadSort['VAR-NT']
  pyReadFilt = pyReadSort[pyReadSort['Pos_Ref_Alt'].isin(shReadSort['Pos_Ref_Alt'])].sort_values(by = ['POSITION']).reset_index().drop(['index'], axis=1)
  #assert pyReadSort.equals(shReadSort)
  if pyReadFilt.equals(shReadSort):
    print("lowfreq final table pass")
  else:
    print("Initial lowfreq final table test FAIL")
    for i in range(len(pyReadFilt.index)):
      print(i)
      pyList = pyReadFilt.loc[i, :].values.flatten().tolist()
      shList = shReadSort.loc[i, :].values.flatten().tolist()
      if pyList != shList:
        print(pyList)
        print(shList)
        print('____________')
    
test_lofreqFinal(df1 = "sample_lofreq-output.tsv", df2 = "test_rose/final-calls/GA-EHC-2884X_L1_lofreq-output.txt")
