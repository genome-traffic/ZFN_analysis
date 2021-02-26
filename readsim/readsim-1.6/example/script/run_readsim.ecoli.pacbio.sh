#!/bin/bash

#set -xv 


cd ../ecoli

for l in 5000; do 
  for c in 5; do 

    ../../src/readsim.py sim fa \
    --ref NC_000913.fna \
    --pre NC_000913.pacbio.reads \
    --rev_strd on \
    --tech pacbio --read_mu $l --cov_mu $c  

  done;
done;


