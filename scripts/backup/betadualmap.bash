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
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.out
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.log
rm -f /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone/*

for f in /media/nikiwind/DATADRIVE2/ZFN_analysis/align/*.fastq
do
	echo " ----------------------------------------------------------------------- Aligning... - $f -------------------------------------------------------------------------------------------------------"

	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa -t 20 -x sensitive -Z -a anchor -d $f -o $f.sam
	grep '^zfn' $f.sam | cut -f1,3,14 > $f.paf
	sed -i -e 's/AS:i://g' $f.paf
	sort -k 3 -nr $f.paf > tmpout && mv tmpout $f.paf

      	refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")
	
	if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
		touch $f.failed
		continue;
	fi

	printf "%s\n" "$refMQ" >> $f.log

	grep "$refvar" /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_start.fa -A 1 > $f.fasta
	for i in {1..8};do cat $f.fasta; done > tmpout && mv tmpout $f.fasta

	for i in `seq 2 16`;
        do	
       		
		if [ $((i%2)) -eq 0 ];
		then
    			echo "-------------------------------------------------------------------------even";
			paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/evenmods.fa > tmpout && mv tmpout $f.fasta
				
			/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 20--rebuild-index --fly-index -x sensitive -Z -a anchor -d $f -o $f.sam
			grep '^zfn' $f.sam | cut -f1,3,14 > $f.paf
			sed -i -e 's/AS:i://g' $f.paf
			sort -k 3 -nr $f.paf > tmpout && mv tmpout $f.paf

      			refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
			refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")


			if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
				head -1 $f.paf > tmpout && mv tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > tmpout && mv tmpout $f.paf
			fi
			
			printf "%s\n" "$refMQ" >> $f.log
			grep "$refvar" $f.fasta -A 1 > tmpout && mv tmpout $f.fasta
			for i in {1..18};do cat $f.fasta; done > tmpout && mv tmpout $f.fasta

		else
    			echo "-------------------------------------------------------------------------odd";
			paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods.fa > tmpout && mv tmpout $f.fasta

			/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 20--rebuild-index --fly-index -x sensitive -Z -a anchor -d $f -o $f.sam
      			grep '^zfn' $f.sam | cut -f1,3,14 > $f.paf
			sed -i -e 's/AS:i://g' $f.paf
			sort -k 3 -nr $f.paf > tmpout && mv tmpout $f.paf

      			refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
			refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")

			if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
				head -1 $f.paf > tmpout && mv tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > tmpout && mv tmpout $f.paf
			fi

			printf "%s\n" "$refMQ" >> $f.log
			grep "$refvar" $f.fasta -A 1 > tmpout && mv tmpout $f.fasta
			for i in {1..16};do cat $f.fasta; done > tmpout && mv tmpout $f.fasta

		fi
		
        done
	
	paste --delimiter='' $f.fasta /media/nikiwind/DATADRIVE2/ZFN_analysis/refs/oddmods_end.fa > tmpout && mv tmpout $f.fasta

	/media/nikiwind/SSD2/bio_software/graphmap/bin/Linux-x64/graphmap align -r $f.fasta -t 20--rebuild-index --fly-index -x sensitive -Z -a anchor -d $f -o $f.sam

	grep '^zfn' $f.sam | cut -f1,3,14 > $f.paf
	sed -i -e 's/AS:i://g' $f.paf
	sort -k 3 -nr $f.paf > tmpout && mv tmpout $f.paf

      	refvar=$(awk -v line="1" -v col="2" 'NR == line { print $col }' < "$f.paf")
	refMQ=$(awk -v line="1" -v col="3" 'NR == line { print $col }' < "$f.paf")

	if ! [[ "$refMQ" =~ ^[0-9]+$ ]]; then
				echo ">>>> switching to minimap";		
				/media/nikiwind/SSD2/bio_software/minimap2/minimap2 -x map-ont $f.fasta $f > $f.paf --secondary=no --for-only
      				sort -k 11 -nr $f.paf > tmpout && mv tmpout $f.paf
				head -1 $f.paf > tmpout && mv tmpout $f.paf
				
				refMQ=$(awk -v line="1" -v col="12" 'NR == line { print $col }' < "$f.paf")
				refvar=$(awk -v line="1" -v col="6" 'NR == line { print $col }' < "$f.paf")
				cut -f6,12 $f.paf > tmpout && mv tmpout $f.paf
	fi


	printf "%s\n" "$refMQ" >> $f.log

	grep "$refvar" $f.fasta -A 1 > tmpout && mv tmpout $f.fasta

	echo "Done aligning... - $f"
	mv $f.* /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
	done

cd /media/nikiwind/DATADRIVE2/ZFN_analysis/aligndone
grep -H ">" *.fasta | column -t > /media/nikiwind/DATADRIVE2/ZFN_analysis/output_fwd.txt


