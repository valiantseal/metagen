conda activate bactopia

cd bactopia_gtdbtk/klebsiella-pneumoniae

mkdir klebOut

bactopia --wf kleborate \
--bactopia ./output \
--outdir ./klebOut \


aws s3 cp --recursive ./klebOut \
s3://transfer-files-emory/ICMC/Lucy/cre_isolates/custom_output/pangenome/klebsiella-pneumoniae/kleborate_output/