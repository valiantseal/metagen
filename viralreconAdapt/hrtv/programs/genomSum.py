import subprocess
import pandas as pd
import glob
import os

workflow = "HRTV"

os.makedirs('summary', exist_ok=True)
os.makedirs('all_summaries', exist_ok=True)
os.chdir('./summary')

samplesDf = pd.read_csv('../input.csv')
samples = samplesDf['sample'].to_list()

def getHrtv():
  dirsList = glob.glob('../output_*')
  for i in dirsList:
    for sample in samples:
      dirName = i.replace('../output_', '')
      inFile = i + '/variants/bowtie2/' + sample + '.ivar_trim.sorted.bam'
      if os.path.isfile(inFile) == True:
        mos_str = 'mosdepth ' +  dirName + '_' + sample + ' ' + inFile
        try:
          subprocess.run(mos_str , shell=True)
        except Exception as e:
          print(e)
        cov_str = 'samtools coverage ' +  inFile + ' > ' + dirName + '_' + sample  + '.coverage.tsv'
        try:
          subprocess.run(cov_str , shell=True)
        except Exception as e:
          print(e)

def getDengue():
  for sample in samples:
    inFile = '../output/variants/bowtie2/' + sample + '.ivar_trim.sorted.bam'
    mos_str = 'mosdepth ' + sample + ' ' + inFile
    subprocess.run(mos_str , shell=True)
    cov_str = 'samtools coverage ' +  inFile + ' > ' + sample  + '.coverage.tsv'
    subprocess.run(cov_str , shell=True)

if workflow == "Dengue":
  getDengue()
elif workflow == "HRTV":
  getHrtv()

os.chdir('../')

# get depths summary
def combMos():
  mosList = glob.glob('./summary/*summary.txt')
  mosSum = pd.DataFrame()
  for i in mosList:
    try:
      mos = pd.read_table(i)
      curFile = i.replace("./summary/", "").replace(".mosdepth.summary.txt", "")
      mos["File"] = curFile
      mosFilt = mos[mos["chrom"] == "total"]
      mosSum = pd.concat([mosSum, mosFilt], axis = 0)
    except Exception as e:
      print
  return(mosSum)

mosSum = combMos()

# get coverage summary
def combCov():
  covList = glob.glob('./summary/*.coverage.tsv')
  covSum = pd.DataFrame()
  for i in covList:
    try:
      cov = pd.read_table(i)
      curFile = i.replace("./summary/", "").replace(".coverage.tsv", "")
      cov["File"] = curFile
      covSum = pd.concat([covSum, cov], axis = 0)
    except Exception as e:
      print(e)
  covSum.columns.values[0] = "ref_name"
  return(covSum)

covSum = combCov()
  
# save 
# run Name
fullDir = os.getcwd()
dirList = fullDir.split('/')
run_name = dirList[6]

mosOut = "all_summaries/depthSummary_" + run_name + ".csv"
mosSum.to_csv(path_or_buf = mosOut, index = False)

covOut = "all_summaries/coverageSummary_" + run_name + ".csv"
covSum.to_csv(path_or_buf = covOut, index = False)
