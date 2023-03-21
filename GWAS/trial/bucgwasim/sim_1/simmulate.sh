conda activate BacGWASim

BacGWASim --num-species 2000 \
--genome-length 30000 \
--mutation-rate 0.06 \
--recomb-rate 0.01 \
--maf 0.01 \
--phen-type "quant" \
--num-causal-var 16 \
--causal-ld-max 0.5 \
--phen-replication 5 \
--cores 14

