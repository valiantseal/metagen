
python ./programs/CltTest.py --help

python ./programs/CltTest.py -i 'GA-EHC-2884X_L1_bbmap-1_lofreq-output.txt' -r 'GA-EHC-2884X_L1_bbmap-1_pos-filter (1).txt' `
-c 0.5 -a 0.2 

python programs/filterLofreqCLT.py -i 'GA-EHC-2884X_L1_bbmap-1_lofreq-output.txt' -r 'GA-EHC-2884X_L1_bbmap-1_pos-filter (1).txt' `
-c 0.5 -a 0.2 

python programs/testResults.py