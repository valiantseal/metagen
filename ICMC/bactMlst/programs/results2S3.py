import subprocess

cmd_str = "aws s3 cp --recursive ./custom_output/ s3://transfer-files-emory/ICMC/Sarah/hflu_2023-03-08/custom_output/"

subprocess.run(cmd_str, shell=True)

