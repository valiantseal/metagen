import subprocess
import pkg_resources

required = {'pandas', 'argparse', 'multiprocessing', 'time', 'math', 'functools'}
installed = {pkg.key for pkg in pkg_resources.working_set}
missing = required - installed

# install required packages
if missing:
    python = sys.executable
    subprocess.check_call([python, '-m', 'pip', 'install', *missing], stdout=subprocess.DEVNULL)
    
