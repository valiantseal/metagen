
sh ./programs/prepInput.sh

for i in $(cat newdir.list); do cd trimKrUnVipr/"$i"; bash -i ../../programs/trimMerge.sh; cd ../../; done

for i in $(cat newdir.list); do cd trimKrUnVipr/"$i"; bash -i ../../programs/FqToFa.sh; cd ../../; done

for i in $(cat newdir.list); do cd trimKrUnVipr/"$i"; bash -i ../../programs/krakenUniq.sh; cd ../../; done

for i in $(cat newdir.list); do cd trimKrUnVipr/"$i"; bash -i ../../programs/blastViprKrakUniq.sh; cd ../../; done
for i in $(cat newdir.list); do cd trimKrUnVipr/"$i"; bash -i ../../programs/selVirus.sh; cd ../../; done


# no filtering kraken uniq vipr
sh ./programs/prepInputNoFiltr.sh
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/merge.sh; cd ../../; done
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/FqToFa.sh; cd ../../; done
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/krakenUniq.sh; cd ../../; done
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/blastViprKrakUniq.sh; cd ../../; done
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/selVirus.sh; cd ../../; done
# try blast without maximum target seq
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/blastViprKrakUniq.sh; cd ../../; done
for i in $(cat newdir.list); do cd krUnVipr/"$i"; bash -i ../../programs/selVirus.sh; cd ../../; done


# check results after fastp filtering
for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../programs/selVirus.sh; cd ../../; done
for i in $(cat newdir.list); do cd work/"$i"; bash -i ../../programs/selVirusKrak2.sh; cd ../../; done

# experiment with splitting sequences:
cd /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads/
conda activate seqkit
seqkit split krakUniq_classified.reads -s 10000 -j 8 -O ./splitSeq10K/
cd splitSeq10K
for i in *.reads; do mkdir "$i"_dir; mv "$i" ./"$i"_dir/reads.fas; done 


cd /home/ubuntu/extraVol/blast_db/nt_data/data_v4

time blastn -db nt -query /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads/splitSeq10K/krakUniq_classified.part_001.reads_dir/reads.fas \
-num_threads 6 -out /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads/splitSeq10K/krakUniq_classified.part_001.reads_dir/NtV4_blast.results \
-outfmt '7 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-max_target_seqs 10

cd /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads
seqkit split krakUniq_classified.reads -s 2000 -j 8 -O ./splitSeq2K/
cd splitSeq2K
for i in *.reads; do mkdir "$i"_dir; mv "$i" ./"$i"_dir/reads.fas; done 

cd /home/ubuntu/extraVol/blast_db/nt_data/data_v4
time blastn -db nt -query /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads/splitSeq2K/krakUniq_classified.part_001.reads_dir/reads.fas \
-num_threads 6 -out /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001/merged_reads/splitSeq2K/krakUniq_classified.part_001.reads_dir/NtV4_blast.results \
-outfmt '7 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
-max_target_seqs 10

# workflow for parallel processing
cd /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-002A_S1_L001
mkdir all_classified_reads
for i in */krakUniq_classified.reads; do cat "$i" >> all_classified_reads/all_classified_reads.fas; done
wc -l */krakUniq_classified.reads
cd all_classified_reads; wc -l all_classified_reads.fas
conda activate seqkit
seqkit stats all_classified_reads.fas
mkdir splitSeq10K
seqkit split all_classified_reads.fas -s 10000 -j 8 -O ./splitSeq10K/
cd splitSeq10K
for i in *.fas; do mkdir "$i"_dir; mv "$i" ./"$i"_dir/reads.fas; done
date >> date.txt; time ls -d */ | parallel -j 3 'cd {} && sh ../../../../../programs/blastNtV4Par.sh'; date >> date.txt # 9 hours 16 CPUs
# 003 sample 98 CPUs # 28 processes completed sucsessfully 58 minutes almost run out of memory cpu usage was adequate
cd /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-003A_S2_L001
mkdir all_classified_reads
for i in */krakUniq_classified.reads; do echo "$i"; cat "$i" >> all_classified_reads/all_classified_reads.fas; done
wc -l */krakUniq_classified.reads
cd all_classified_reads; wc -l all_classified_reads.fas
conda activate seqkit
seqkit stats all_classified_reads.fas
mkdir splitSeq10K
seqkit split all_classified_reads.fas -s 10000 -j 20 -O ./splitSeq10K/
cd splitSeq10K
for i in *.fas; do mkdir "$i"_dir; mv "$i" ./"$i"_dir/reads.fas; done
time ls -d */ | parallel -j 28 'cd {} && sh ../../../../../programs/blastNtV4Par.sh'

# sample B2E22-004A_S3_L001  98 CPUs 28 processes completed in 84 minutes
cd /home/ubuntu/extraVol/metagenClass/2022-02-25/krUnVipr/B2E22-004A_S3_L001
mkdir all_classified_reads
for i in */krakUniq_classified.reads; do echo "$i"; cat "$i" >> all_classified_reads/all_classified_reads.fas; done
wc -l */krakUniq_classified.reads
cd all_classified_reads; wc -l all_classified_reads.fas
conda activate seqkit
seqkit stats all_classified_reads.fas
mkdir splitSeq10K
seqkit split all_classified_reads.fas -s 10000 -j 20 -O ./splitSeq10K/
cd splitSeq10K
for i in *.fas; do mkdir "$i"_dir; mv "$i" ./"$i"_dir/reads.fas; done
time ls -d */ | parallel -j 28 'cd {} && sh ../../../../../programs/blastNtV4Par.sh'

# sum blast Nt output
sh ./programs/sumVirusBlastNt.sh