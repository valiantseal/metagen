bash -i ./programs/prepInput.sh

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/filterMerge.sh; cd ../../; done

for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../../programs/faKraken2.sh; cd ../../; done