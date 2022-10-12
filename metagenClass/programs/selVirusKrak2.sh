cat sample.name >> ../selectVirusKrak2.summary

for j in $(cat merged.list)
do 
cd "$j"
mkdir select_viruses
rm select_viruses/Report_alphaherpesvirus.krak2
rm select_viruses/Report_mastadenovirus.krak2
rm select_viruses/Report_polyomavirus.krak2

rm select_viruses/ViprBlast_results_alphaherpesvirus.krak2
rm select_viruses/ViprBlast_results_mastadenovirus.krak2
rm select_viruses/ViprBlast_results_polyomavirus.krak2

grep -i "alphaherpesvirus 2" sample.report > select_viruses/Report_alphaherpesvirus.krak2
grep -i "Mastadenovirus D" sample.report > select_viruses/Report_mastadenovirus.krak2
grep -i "JC Polyomavirus" sample.report > select_viruses/Report_polyomavirus.krak2

grep -i "alphaherpesvirus 2" blast.results > select_viruses/ViprBlast_results_alphaherpesvirus.krak2
grep -i "Mastadenovirus D" blast.results > select_viruses/ViprBlast_results_mastadenovirus.krak2
grep -i "JC Polyomavirus" blast.results > select_viruses/ViprBlast_results_polyomavirus.krak2

cd select_viruses
echo "$j" >> ../../../selectVirusKrak2.summary

for file in *.krak2
do
wc -l "$file" >> ../../../selectVirusKrak2.summary
done
cd ../

cd ../
done