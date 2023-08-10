import os
import fnmatch
import pandas as pd
from multiprocessing import Pool, freeze_support

os.chdir("C:/Users/abomb/OneDrive - Emory University/virus/iSNVs")

def getSamplesNames(df, curType):
  allSampleNames = []
  curDf = pd.read_csv(df)
  for i in range(len(curDf.index)):
    if curType == "metaseq":
      curSeq = curDf.loc[i, "MetaseqNames"]
    elif curType == "ampseq":
      curSeq = curDf.loc[i, "Missing_Ampseq"]
    if curSeq  != "NA" and not pd.isna(curSeq):
      curSample = curDf.loc[i, "Sample_id"]
      allSampleNames.append(curSample)
  return allSampleNames

def findFiles(curSample):
    combPaths = []
    combSizes = []
    combNames = []
    root_directory = "Z:/piantadosilab/EICC_Output/viralrecon/EHC_C19_ampseq/"
    ignored_directories = ["work", "water_samples", "results", "pfizer_tsv", "mosdepth_report_data"]
    pattern = "*" + curSample.replace("-", "*").replace("_", "*") + "*"
    print(pattern)
    for root, dirs, files in os.walk(root_directory):
        # Remove ignored directory names from the list of directories
        dirs[:] = [d for d in dirs if d not in ignored_directories]
        
        matching_files = fnmatch.filter(files, pattern)
        for matching_file in matching_files:
            file_path = os.path.join(root, matching_file)
            file_size = round(os.path.getsize(file_path) / (1024 * 1024), 2)
            combPaths.append(file_path)
            combSizes.append(file_size)
            combNames.append(curSample)
    if len(combPaths) > 0 and len(combSizes) > 0:
      try:
        combDat = pd.DataFrame({"Sample_id":combNames, "Pbnas_path":combPaths, "File_size":combSizes})
        outFile = f"pbnas_ampseq_paths/{curSample}.csv"
        combDat.to_csv(path_or_buf=outFile, index=False)   
      except:
         pass

if __name__ == '__main__':
    curSamples = getSamplesNames(df="SeqSamples_Eval_Depth10.csv", curType="ampseq")

    # Add this line to support multiprocessing on Windows
    freeze_support()

    with Pool(10) as pool:
        pool.map(findFiles, curSamples)


#findFiles(curSample = "EHC-C19-3838P")