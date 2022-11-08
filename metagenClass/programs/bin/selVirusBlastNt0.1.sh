rm *.par

cat ../../../../programs/virus.list | parallel -j 27 'grep -i {} NtV4_blast.results > {}.par'