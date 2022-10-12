#grep -i "alphaherpesvirus 2" NtV4_blast.results > alphaherpesvirus.sel
#grep -i "Mastadenovirus" NtV4_blast.results > mastadenovirus.sel
#grep -i "Polyomavirus" NtV4_blast.results > polyomavirus.sel

#grep -i "herpesvirus 2" NtV4_blast.results > herpesvirus2.sel

#grep -i "herpes simplex virus 2" NtV4_blast.results > herpessimplex2.sel

rm *.sel

for i in $(cat ../../../../virus.list)
do
echo "$i" > current.vir
sed -i "s/ /_/g" current.vir
virus=$(cat current.vir)
grep -i "$i" NtV4_blast.results > "$virus".sel
done
