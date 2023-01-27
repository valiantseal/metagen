conda activate seqtk
cd testReads

seqtk subseq /home/ubuntu/extraVol/metagenClass/Nima_p2/process/Water-042022_S14_L001/krakUniq_classified.reads \
waterReads.list > waterReads.fasta

grep -F -f waterReads.list \
/home/ubuntu/extraVol/metagenClass/Nima_p2/blastNtSummary/target_results/Water-042022_S14_L001__Human_alphaherpesvirus_1 \
> waterMatchBlast.res 

grep -i "Human alphaherpesvirus 1" /home/ubuntu/extraVol/metagenClass/Nima_p2/process/Water-042022_S14_L001/krakUniq_sample.report \
> waterKraken.report
