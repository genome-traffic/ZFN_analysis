#!/bin/bash

MQcutoff=61

#first cut in pieces
#count reads
#awk '{s++}END{print s/4}' file.fastq
#split
#perl fastq-splitter.pl f0.fastq --n-parts 4001 --measure count

#https://lh3.github.io/minimap2/minimap2.html

rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.sam
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.paf
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.ref
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone/*

for f in /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.fastq
do
	echo " ----------------------------------------------------------------------- Aligning... - $f -------------------------------------------------------------------------------------------------------"

	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa $f > $f.paf --secondary=no --for-only
	
	sort -k 12 -n $f.paf > tmpout && mv tmpout $f.paf
	head -1 $f.paf > tmpout && mv tmpout $f.paf
	
	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
	
	if [ -z "$refMQ" ]; then
	rm $f.paf
	rm $f
	continue
	fi
	
	if [ $refMQ -gt $MQcutoff ]; then
	sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
	fi		
	
	refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
	refL=$(awk -v line="1" -v col="11" 'NR == line { print $col }' < "$f.paf")

	printf "%s\t" "$refL" > $f.log
	printf "%s\n" "$refMQ" >> $f.log

	grep "$refvar" /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa -A 1 > $f.ref
	for i in {1..8};do cat $f.ref; done > tmpout && mv tmpout $f.ref

	for i in `seq 2 16`;
        do	
       		echo $i
		if [ $((i%2)) -eq 0 ];
		then
    		echo "-------------------------------------------------------------------------even";
		paste --delimiter='' $f.ref /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/evenmods.fa > tmpout && mv tmpout $f.ref
		
		/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
		sort -k 12 -n $f.paf > tmpout && mv tmpout $f.paf
		refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
		
		if [ $refMQ -gt $MQcutoff ]; then
      		sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
		fi		

		head -1 $f.paf > tmpout && mv tmpout $f.paf

		refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
		refL=$(awk -v line="1" -v col="11" 'NR == line { print $col }' < "$f.paf")
		
		printf "%s\t" "$refL" >> $f.log
		printf "%s\n" "$refMQ" >> $f.log

		refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
		grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref
		
		for i in {1..18};do cat $f.ref; done > tmpout && mv tmpout $f.ref

		else
    		echo "-------------------------------------------------------------------------odd";
		paste --delimiter='' $f.ref /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods.fa > tmpout && mv tmpout $f.ref  
	
		/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
		sort -k 12 -n $f.paf > tmpout && mv tmpout $f.paf
		refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
		
		if [ $refMQ -gt $MQcutoff ]; then
      		sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
		fi		

		head -1 $f.paf > tmpout && mv tmpout $f.paf

		refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
		refL=$(awk -v line="1" -v col="11" 'NR == line { print $col }' < "$f.paf")
		
		printf "%s\t" "$refL" >> $f.log
		printf "%s\n" "$refMQ" >> $f.log
		
		refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
		grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref

		for i in {1..16};do cat $f.ref; done > tmpout && mv tmpout $f.ref

		fi

        done
	
	paste --delimiter='' $f.ref /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_end.fa > tmpout && mv tmpout $f.ref

	#final alignement	

	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
	sort -k 12 -n $f.paf > tmpout && mv tmpout $f.paf
	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
		
	if [ $refMQ -gt $MQcutoff ]; then
      	sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
	fi		

	head -1 $f.paf > tmpout && mv tmpout $f.paf

	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
	refL=$(awk -v line="1" -v col="11" 'NR == line { print $col }' < "$f.paf")
	
	printf "%s\t" "$refL" >> $f.log
	printf "%s\n" "$refMQ" >> $f.log

	refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
	grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref
		
	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -ax map-ont $f.ref $f > $f.sam  --secondary=no --for-only
	mv $f.ref /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	mv $f.sam /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	mv $f.paf /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	mv $f.log /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	echo "Done aligning... - $f"

done

cd /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
grep -H "f_" *.paf | cut -f1,2,6,10,11,12 | column -t > /media/nikiwind/DATADRIVE2/ZFN_analysis/output_fwd.txt


