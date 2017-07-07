#!/bin/bash
for i in {1..30};
do

mkdir sample${i}-quant

kallisto quant -t 12 -i Homo_sapiens.GRCh38.88.cdna.kallisto.idx -o sample${i}-quant ../data/sample_${i}/sample_${i}_1.fasta ../data/sample_${i}/sample_${i}_2.fasta

done
