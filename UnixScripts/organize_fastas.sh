#!/bin/bash

for i in {1..30};
do
mkdir sample_${i}
mv sample_${i}_1.fasta.gz sample_${i}
mv sample_${i}_2.fasta.gz sample_${i}
done
