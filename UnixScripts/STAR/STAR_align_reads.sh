#!/bin/bash

for i in {13..30};
do
cd sample_${i}
STAR --genomeDir ../../star-index/ --readFilesIn sample_${i}_1.fasta sample_${i}_2.fasta --runThreadN 12 --quantMode TranscriptomeSAM
cd ..
done
