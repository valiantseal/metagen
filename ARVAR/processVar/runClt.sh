conda activate arvar

python ./programs/getVarFilesCLT.py --help

python ./programs/getVarFilesCLT.py -i input/output.sam \
-r references/MN908947.3.fna

time python ./programs/getVarFilesCLT.py -i test_rose/GA-EHC-2884X_L1_hst2.sam \
-r test_rose/SARSCov2_Ref.fasta -t 8 -a 0.001

python ./programs/test_GetVarFiles.py