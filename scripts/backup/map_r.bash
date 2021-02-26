#!/bin/bash


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

	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_endstart.fa $f > $f.paf --secondary=no --for-only
	
	awk 'NR == 1 || $12 < min {line = $0; min = $12}END{print line}' $f.paf > tmpout && mv tmpout $f.paf
	refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
	echo "Best reference is $refvar with a mapping quality of $refMQ !"
	
	if [ -z "$refMQ" ]; then
	rm $f.paf
	#rm $f
	continue

	
	grep "$refvar" /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_endstart.fa -A 1 > $f.ref
	for i in {1..8};do cat $f.ref; done > tmpout && mv tmpout $f.ref
	
	for i in `seq 2 16`;
        do	
       		echo $i
		if [ $((i%2)) -eq 0 ];
		then
    		echo "even";
		paste --delimiter='' /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/evenmods_r.fa $f.ref > tmpout && mv tmpout $f.ref
		sed -i 's/_>/_/g' $f.ref

		/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
		awk 'NR == 1 || $12 < min {line = $0; min = $12}END{print line}' $f.paf > tmpout && mv tmpout $f.paf
		refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
		grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref
		for i in {1..18};do cat $f.ref; done > tmpout && mv tmpout $f.ref

		else
    		echo "odd";
		paste --delimiter='' /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_r.fa $f.ref > tmpout && mv tmpout $f.ref
		sed -i 's/_>/_/g' $f.ref

		/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
		awk 'NR == 1 || $12 < min {line = $0; min = $12}END{print line}' $f.paf > tmpout && mv tmpout $f.paf
		refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
		grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref
		for i in {1..16};do cat $f.ref; done > tmpout && mv tmpout $f.ref

		fi

        done
	
	paste --delimiter='' /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_startend.fa $f.ref > tmpout && mv tmpout $f.ref
	sed -i 's/_>/_/g' $f.ref

	#final alignement	

	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.ref $f > $f.paf --secondary=no --for-only
	awk 'NR == 1 || $12 < min {line = $0; min = $12}END{print line}' $f.paf > tmpout && mv tmpout $f.paf
	refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
	grep "$refvar" $f.ref -A 1 > tmpout && mv tmpout $f.ref
	/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -ax map-ont $f.ref $f > $f.sam  --secondary=no --for-only
	mv $f.ref /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	mv $f.sam /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	mv $f.paf /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone

	echo "Done aligning... - $f"

done

cd /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
grep -H "f_" *.paf | cut -f1,2,6,10,11,12 | column -t > /media/nikiwind/DATADRIVE2/ZFN_analysis/output_rev.txt


