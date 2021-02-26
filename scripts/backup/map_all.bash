#!/bin/bash

for f in /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.fastq
do
	echo " ----------------------------------------------------------------------- Aligning... - $f -------------------------------------------------------------------------------------------------------"

	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/ZFarray_generalized.fa -t 20 --rebuild-index -Z -a anchor -d $f -o $f.sam
	#/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/ZFarray_generalized.fa $f > $f.paf --secondary=no --for-only
	cp /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/ZFarray_generalized.fa $f.fasta
	#echo "Done aligning... - $f"

done

#grep -H "AD" *.paf | cut -f1,6,10,11,12 | column -t > /media/nikiwind/DATADRIVE2/ZFN_analysis/output.txt


