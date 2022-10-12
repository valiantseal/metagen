rm target.viruses

cat ../../../../virus.list | parallel -j 26 'grep -i {} NtV4_blast.results > {}.par'