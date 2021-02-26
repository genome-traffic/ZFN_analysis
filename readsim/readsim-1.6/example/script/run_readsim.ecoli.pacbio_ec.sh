#!/bin/bash

#set -xv 


cd ../ecoli

for l in 5000; do 
  for c in 5; do 

    ../../src/readsim.py sim fa \
    --ref NC_000913.100k.fna \
    --pre NC_000913.100k.pacbio_ec.reads \
    --rev_strd on \
    --tech pacbio_ec --read_mu $l --cov_mu $c --read_dist ../pacbio.p5-c3.len.dist 

  done;
done;


