conda activate miniwdl

cd /home/ubuntu/viral-ngs/myPipe/patient288/metaSeq

Rscript --vanilla ./programs/inputFile.R

for i in $(cat samples.txt); do mkdir -p process_files/"$i"; mv ./input/"$i"* ./process_files/"$i";
echo "$i" >process_files/"$i"/sample.txt;
done

cd process_files
# convert fastq to bam. takes 2 cpus 3G ram per sample reserved 
for i in $(cat ../samples.txt); do cd "$i"; sh ../../programs/uBam.sh; cd ../; done

# assemble genome denovo with filtering 
for i in $(cat ../samples.txt); do cd "$i"; sh ../../programs/assembleDenovo.sh; cd ../; done


# prepare files for manula iSNVs call
for i in $(cat ../samples.txt); do mkdir -p "$i"/iSnvs; cp "$i"/sample.txt "$i"/iSnvs/; done

# filter reads with fastP or trimgalore
for i in $(cat ../samples.txt); do cd "$i"; bash -i ../../programs/trimgalore.sh; cd ../; done


# convert to bam
for i in $(cat ../samples.txt); do cd "$i"/iSnvs; sh ../../../programs/uBam.sh; cd ../../; done

# clean with bmtagger
# human reads 8 cpus (max at the moment) 14GB ram per samples reserved, appears to use less than half
for i in $(cat ../samples.txt); do cd "$i"/iSnvs; sh ../../../programs/bmtaggerHum.sh; cd ../../; done

# allign
for i in $(cat ../samples.txt); do cd "$i"/iSnvs; sh ../../../programs/allign.sh; cd ../../; done

# call iSnv
for i in $(cat ../samples.txt); do cd "$i"/iSnvs; sh ../../../programs/isnv.sh; cd ../../; done


