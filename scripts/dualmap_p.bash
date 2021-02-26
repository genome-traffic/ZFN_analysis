#!/bin/bash

#first cut in pieces
#count reads
#awk '{s++}END{print s/4}' file.fastq
#split
#perl fastq-splitter.pl f0.fastq --n-parts 4001 --measure count
#find . -size 0 -delete



#https://lh3.github.io/minimap2/minimap2.html

find . -name "*fastq.*" -delete

mapitNANO() {
	
	f=${1}
	echo " ----------------------------------------------------------------------- Aligning... - $f -------------------------------------------------------------------------------------------------------"
	
	cd /media/nikiwind/DATADRIVE2/ZFN_analysis/align/
	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa -t 6 --rebuild-index -x sensitive -Z -a anchor -d $f -o $f.sam
	grep 'AS:i:' $f.sam | cut -f1,3,14 > $f.paf
	
	sed -i -e 's/AS:i://g' $f.paf
	sort -k 3 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf
	
	refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")
	
	if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
		
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				head -1 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")

				if [ -z "$refMQ" ]; then
				touch $f.failed
				refMQ="BAD"
				return 1
				fi
				
				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				printf "%s\t" "mmap" >> $f.log

	else
				printf "%s\t" "gmap" >> $f.log
	fi

	printf "%s\n" "$refMQ" >> $f.log

	grep "$refvar" /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa -A 1 > $f.fasta
	for i in {1..16};do cat $f.fasta; done > $f.tmpout && mv $f.tmpout $f.fasta

	
	for i in `seq 2 16`;
        do	
       		
		if [ $((i%2)) -eq 0 ];
		then
    		echo "-------------------------------------------------------------------------even";
		paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/evenmods.fa > $f.tmpout && mv $f.tmpout $f.fasta
				
		/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 6 --rebuild-index -x sensitive -Z -a anchor -d $f -o $f.sam
		grep 'AS:i:' $f.sam | cut -f1,3,14 > $f.paf
		sed -i -e 's/AS:i://g' $f.paf
		sort -k 3 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf

      		refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
		refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")
		
		if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
		
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				head -1 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")

				if [ -z "$refMQ" ]; then
				touch $f.failed
				refMQ="BAD"
				return 1
				fi

				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				printf "%s\t" "mmap" >> $f.log

		else
				printf "%s\t" "gmap" >> $f.log
		fi

		printf "%s\n" "$refMQ" >> $f.log
		grep "$refvar" $f.fasta -A 1 > $f.tmpout && mv $f.tmpout $f.fasta
		for i in {1..18};do cat $f.fasta; done > $f.tmpout && mv $f.tmpout $f.fasta

		else
    		echo "-------------------------------------------------------------------------odd";
		paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods.fa > $f.tmpout && mv $f.tmpout $f.fasta

		/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 6 --rebuild-index -x sensitive -Z -a anchor -d $f -o $f.sam
      		grep 'AS:i:' $f.sam | cut -f1,3,14 > $f.paf
		sed -i -e 's/AS:i://g' $f.paf
		sort -k 3 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf

      		refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
		refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")

		if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
		
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				head -1 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")

				if [ -z "$refMQ" ]; then
				touch $f.failed
				refMQ="BAD"
				return 1
				fi

				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				printf "%s\t" "mmap" >> $f.log
		else
				printf "%s\t" "gmap" >> $f.log
		fi

		printf "%s\n" "$refMQ" >> $f.log
		grep "$refvar" $f.fasta -A 1 > $f.tmpout && mv $f.tmpout $f.fasta
		for i in {1..16};do cat $f.fasta; done > $f.tmpout && mv $f.tmpout $f.fasta

		fi
		
        done
	
	paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_end.fa > $f.tmpout && mv $f.tmpout $f.fasta
	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 6 --rebuild-index -x sensitive -Z -a anchor -d $f -o $f.sam
	grep 'AS:i:' $f.sam | cut -f1,3,14 > $f.paf
	sed -i -e 's/AS:i://g' $f.paf
	sort -k 3 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf

	refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")

	if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
		
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				head -1 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")

				if [ -z "$refMQ" ]; then
				touch $f.failed
				refMQ="BAD"
				return 1
				fi

				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > $f.tmpout && mv $f.tmpout $f.paf
				printf "%s\t" "mmap" >> $f.log
	else
				printf "%s\t" "gmap" >> $f.log
	fi
	
	printf "%s\n" "$refMQ" >> $f.log
	grep "$refvar" $f.fasta -A 1 > $f.tmpout && mv $f.tmpout $f.fasta
	
	rm -f $f.paf
	rm -f $f.fasta.gmidx
	
	echo "Done aligning... - $f"



	}
export -f mapitNANO

ls /media/nikiwind/DATADRIVE2/ZFN_analysis/align/ | grep \.fastq | parallel "mapitNANO /media/nikiwind/DATADRIVE2/ZFN_analysis/align/{}"

cd /media/nikiwind/DATADRIVE2/ZFN_analysis/align
grep -H ">" *.fasta | column -t > /media/nikiwind/DATADRIVE2/ZFN_analysis/output_dualmap.txt

#gedit format into tab seperated columns
#awk 'NF==18' output_dualmap_run2.txt > run2.txt


