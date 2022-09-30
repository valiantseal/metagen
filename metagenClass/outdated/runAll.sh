bash -i ./programs/prepInput.sh

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/filterMerge.sh; cd ../../; done

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/FqToFa.sh; cd ../../; done

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/kraken2.sh; cd ../../; done

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/prepForBlast.sh; cd ../../; done

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/blastVipr.sh; cd ../../; done