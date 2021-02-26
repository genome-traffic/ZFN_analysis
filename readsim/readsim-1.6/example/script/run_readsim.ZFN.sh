#!/bin/bash

#set -xv 


cd ../ZFN

for l in 6200; do 
  for c in 20; do 

	#../../src/readsim.py sim fq --ref ZFN.fna --pre zfn.nanopore.reads --rev_strd off --tech nanopore --read_mu $l --cov_mu $c --read_dist normal 
	../../src/readsim.py sim fq --ref ZFN.fna --pre zfn.nanopore.reads --rev_strd off --tech nanopore --read_mu $l --cov_mu $c --read_dist normal
	#../../src/readsim.py sim fq --ref ZFNrev.fna --pre zfn.nanopore.reads --rev_strd off --tech perfect --read_mu $l --cov_mu $c --read_dist normal
	#../../src/readsim.py sim fq --ref ZFNtrunc.fna --pre zfn.nanopore.reads --rev_strd off --tech perfect --read_mu $l --cov_mu $c --read_dist normal 
	#../../src/readsim.py sim fq --ref random.fna --pre zfn.nanopore.reads --rev_strd off --tech nanopore --read_mu $l --cov_mu $c --read_dist normal
  done;
done;


