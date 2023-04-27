from Bio import SeqIO
import subprocess
import shutil
import os

      
def editFasta(inFile, outFile):
  with open(inFile, "r") as in_handle, open(outFile, "w") as out_handle:
    for record in SeqIO.parse(in_handle, "fastq"):
      record.id = record.id.replace("/1", "")
      record.id = record.id.replace("/2", "")
      record.description = record.description.replace("/1", "")
      record.description = record.description.replace("/2", "")
      SeqIO.write(record, out_handle, "fastq")


editFasta(inFile = "Short_read_test_HRTV_L1.fastq", outFile = "Short_read_test_HRTV_L_R1.fastq")    

editFasta(inFile = "Short_read_test_HRTV_L2.fastq", outFile = "Short_read_test_HRTV_L_R2.fastq") 

os.makedirs("input", exist_ok = True)

def zipFile(infile):
  cmd_str = f'gzip {infile}'
  subprocess.run(cmd_str, shell = True)

zipFile(infile = "Short_read_test_HRTV_L_R1.fastq")
zipFile(infile = "Short_read_test_HRTV_L_R2.fastq")

shutil.copy("Short_read_test_HRTV_L_R2.fastq.gz", "input/")
shutil.copy("Short_read_test_HRTV_L_R1.fastq.gz", "input/")
