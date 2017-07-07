#!/bin/bash
for i in {1..30};
do
mkdir sample_${i}
cd sample_${i}
wget http://rafalab.rc.fas.harvard.edu/rnaseqcomp/simulation/ORIGINAL_READS/sample_${i}_1.fasta.gz
wget http://rafalab.rc.fas.harvard.edu/rnaseqcomp/simulation/ORIGINAL_READS/sample_${i}_2.fasta.gz
cd ..
done
