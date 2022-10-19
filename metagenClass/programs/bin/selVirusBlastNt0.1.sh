rm *.par

cat ../../../../virus.list | parallel -j 27 'grep -i {} NtV4_blast.results > {}.par'