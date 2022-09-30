cat sample.name >> ../selectVirus.summary

for j in $(cat merged.list)
do 
cd "$j"
mkdir select_viruses
rm select_viruses/krakUniqReport_alphaherpesvirus.sel
rm select_viruses/krakUniqReport_mastadenovirus.sel
rm select_viruses/krakUniqReport_polyomavirus.sel

rm select_viruses/ViprBlast_results_alphaherpesvirus.sel
rm select_viruses/ViprBlast_results_mastadenovirus.sel
rm select_viruses/ViprBlast_results_polyomavirus.sel

grep -i "alphaherpesvirus 2" krakUniq_sample.report > select_viruses/krakUniqReport_alphaherpesvirus.sel
grep -i "Mastadenovirus D" krakUniq_sample.report > select_viruses/krakUniqReport_mastadenovirus.sel
grep -i "JC Polyomavirus" krakUniq_sample.report > select_viruses/krakUniqReport_polyomavirus.sel

grep -i "alphaherpesvirus 2" krakenUniq_Vipr_blast.results > select_viruses/ViprBlast_results_alphaherpesvirus.sel
grep -i "Mastadenovirus" krakenUniq_Vipr_blast.results > select_viruses/ViprBlast_results_mastadenovirus.sel
grep -i "Polyomavirus" krakenUniq_Vipr_blast.results > select_viruses/ViprBlast_results_polyomavirus.sel

cd select_viruses
echo "$j" >> ../../../selectVirus.summary

for file in *.sel
do
wc -l "$file" >> ../../../selectVirus.summary
done
cd ../

cd ../
done
