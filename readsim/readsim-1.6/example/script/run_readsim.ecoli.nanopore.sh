#!/bin/bash

#set -xv 


cd ../ecoli

for l in 50000; do 
  for c in 5; do 

    ../../src/readsim.py sim fa \
    --ref NC_000913.fna \
    --pre NC_000913.nanopore.reads \
    --rev_strd on \
    --tech nanopore --read_mu $l --cov_mu $c --read_dist normal 

  done;
done;


