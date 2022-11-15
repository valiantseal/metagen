rm *.par

cat ../../../../programs/virus.list | parallel -j 51 'grep -i {} NtV4_blast.results > {}.par'