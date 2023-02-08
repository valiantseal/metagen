conda activate seqtk


cd process; ls -d */ | parallel -j 4 'cd {} && seqtk subseq ./krakUniq_classified.reads \
krakenReads.id > selectKraken.reads'; cd ../