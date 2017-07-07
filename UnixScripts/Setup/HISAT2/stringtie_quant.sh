#!/bin/bash

for i in {2..30};
do
stringtie -p 8 -G gtf/Homo_sapiens.GRCh38.88.gtf -e -B -o sample_${i}/sample_${i}_transcripts.gtf -A sample_${i}/sample_${i}_gene_abundances.tsv sample_${i}/sample_${i}_alns.sorted.bam
done

