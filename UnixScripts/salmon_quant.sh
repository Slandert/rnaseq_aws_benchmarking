#!/bin/bash
for fn in data/sample_{1..30};
do
samp=`basename ${fn}`
echo "Processing sample ${samp}"

salmon quant -i Homo_sapiens.GRCh37.75_quasi_index/ -l IU -1 data/${samp}_1.fasta.gz -2 data/${samp}_2.fasta.gz -p 8 -o quants/


done
