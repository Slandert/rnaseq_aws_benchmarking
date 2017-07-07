#!/bin/bash

for i in {21..30};
do

hisat2 -p 8 \
-f \
--dta \
--known-splicesite-infile /home/ubuntu/gtf/splicesites.txt \
-x /home/ubuntu/genome/hs-index/hs-index \
-1 /home/ubuntu/fastas/sample_${i}/sample_${i}_1.fasta.gz \
-2 /home/ubuntu/fastas/sample_${i}/sample_${i}_2.fasta.gz \
-S ./sample_${i}_align.sam
done

