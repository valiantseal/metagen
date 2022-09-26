conda activate bactopia-dev
mkdir -p bactopia_summary/combined

for i in $(cat bacteria.list)
do
cd ./bactopia_gtdbtk/"$i"
bactopia-summary output --outdir ./bactopia_summary
cp ./bactopia_summary/bactopia-report.txt ../../bactopia_summary/"$i".tsv
cd ../../
done

Rscript --vanilla ./programs/combSummary.R




