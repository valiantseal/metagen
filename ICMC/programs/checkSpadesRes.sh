for i in $(cat newdir.txt); do
if [ -f ./gtdbtk/input/"$i".fasta ]; then
    echo "$i exists" >> ./gtdbtk/completeSpades.txt
else 
    echo "$i  failed." >> ./gtdbtk/failedSpades.txt
fi
done

echo 'Successful'
wc -l ./gtdbtk/completeSpades.txt

echo 'Failed'
wc -l  ./gtdbtk/failedSpades.txt
