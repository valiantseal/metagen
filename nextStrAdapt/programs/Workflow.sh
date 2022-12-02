cd /home/ubuntu/strain/ncov
snakemake --profile my_profiles/georgia021522Global -p
# emory all sequences 03/07/22 with north america reference 
sh ./programs/emoryCov.sh

# Ludy's MLS data may 2022
snakemake --profile my_profiles/ludyMlsGlob5422 -p
# things to correct:
# format date
# chnage metadata id to lower case
# remove second name for each sequence in fasta file

# visualize outside emory
nextstrain view auspice --host "10.66.123.250" --port 80

#outside emory VPN
http://170.140.124.37:80/ncov/penCovMayGlob
170.140.124.37


# data for shiny App
snakemake --profile my_profiles/forShiny -p


snakemake --profile my_profiles/example -p

# run for ludy
snakemake --profile my_profiles/ludyMlsLoc063022 -p


