conda activate BacGWASim

BacGWASim --num-species 2000 \
--genome-length 2800000 \
--mutation-rate 0.06 \
--recomb-rate 0.01 \
--maf 0.02 \
--phen-type "quant" \
--num-causal-var 16 \
--causal-maf-min 0.2 \
--causal-ld-max 0.5 \
--phen-replication 10 \
--cores 14

