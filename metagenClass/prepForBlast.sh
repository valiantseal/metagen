for j in $(cat merged.list)
do 
cd "$j"

sample=$(cat merge.type)

awk -F '\t' '{ if  ( $1 == "C" ) { print }}' sample.kraken  > classified.reads


awk -F '\t' '{ print $2 }' classified.reads > classified_reads.names


conda activate seqtk
seqtk subseq "$sample".fa classified_reads.names > readsForBlast.fa
cd ../

done

