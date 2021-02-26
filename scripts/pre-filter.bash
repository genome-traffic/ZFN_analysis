#!/bin/bash

for f in /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.fastq
do
	echo " ----------------------------------------------------------------------- Aligning... - $f -------------------------------------------------------------------------------------------------------"

	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/ZFarray_generalized.fa -t 38 -Z -a anchor -d $f -o $f.sam
	cp /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/ZFarray_generalized.fa $f.fasta
	samtools fastq $f.sam -F 4 > $f.mapping.fastq

done



