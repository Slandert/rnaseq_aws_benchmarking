#!/bin/bash
for i in {1..30};
do
cd sample_${i}

rsem-calculate-expression --bam --no-bam-output -p 12 --no-qualities --paired-end Aligned.toTranscriptome.out.bam ~/rsem_ref/ref ./Quant

cd ..
done
