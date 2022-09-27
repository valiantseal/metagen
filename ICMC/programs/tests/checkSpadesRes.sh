for i in $(cat newdir.txt); do
if [ -f ./gtdbtk/input/"$i".fasta ]; then
    echo "$i exists" >> ./gtdbtk/completeSpades.txt
else 
    echo "$i  failed." >> ./gtdbtk/failedSpades.txt
fi
done

echo ''
echo 'Successful'
wc -l ./gtdbtk/completeSpades.txt
echo ''

echo ''
echo 'Failed'
FILE=./gtdbtk/failedSpades.txt
if [ -f "$FILE" ]; then
    wc -l  ./gtdbtk/failedSpades.txt
else 
    echo 0
fi

