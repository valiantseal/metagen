
parNumb=$(wc -l < ../../lenFiltVir.reads)

#echo "$parNumb"

rm -rf ./virReadsFiltLen

mkdir -p ./virReadsFiltLen

#cat ../../lenFiltVir.reads | parallel -j "$parNumb" 'grep -F {} ./NtV4_blast.results >> ./virReadsFiltLen/{}.par'

cat ../../lenFiltVir.reads | parallel -j 51 'grep -F {} ./NtV4_blast.results >> ./virReadsFiltLen/{}.par'
