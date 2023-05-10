import os
import subprocess

threads = '8'
mapqual = '20'
basequal = '20'
refFasta = './references/MN908947.3.fna'

def lofreq(bam, ref_seq, t):
  cmd_str = f"lofreq call-parallel --pp-threads {t} -f {ref_seq} -o sample_lf.vcf {bam} 2>sample_lofreq-log.txt"
  subprocess.run(cmd_str, shell = True)

lofreq(bam = 'sample_dindel.tmp.bam', ref_seq = refFasta, t = threads)


