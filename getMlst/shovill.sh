conda activate shovill


shovill --outdir ./shovill_out --R1 trimmed_val_1.fq.gz --R2 trimmed_val_2.fq.gz \
--cpus 16 --assembler "spades" --force