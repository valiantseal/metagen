import os
import shutil
import subprocess
import glob

def getProjectName():
    curDir = os.path.basename(os.getcwd())
    return curDir


def downloadProject(projName):
    os.makedirs("project_files", exist_ok=True)
    cmd_str = f'bs download project --name {projName} --output project_files'
    subprocess.run(cmd_str, shell=True)

def moveFiles():
    os.makedirs("fastqs", exist_ok=True)
    file_list = glob.glob("project_files/**/*.fastq.gz", recursive=True)
    destination_directory = "fastqs/"
    waterDir = "water_fastqs"
    os.makedirs(waterDir, exist_ok=True)
    for file_path in file_list:
        if "water" in file_path.lower():
            shutil.move(file_path, waterDir)
        else:
            shutil.move(file_path, destination_directory)

if __name__ == "__main__":
    curProject = getProjectName()
    downloadProject(projName = curProject)
    moveFiles()


