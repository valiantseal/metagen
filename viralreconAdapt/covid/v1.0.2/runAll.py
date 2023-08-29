import subprocess

subprocess.run("programs/download.py", shell=True)

subprocess.run("programs/inputFile.py", shell=True)

subprocess.run("programs/runViralrecon.py", shell=True)

subprocess.run("programs/results2S3.py", shell=True)