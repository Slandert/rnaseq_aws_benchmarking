#!/bin/bash

for i in {21..30};
do
samtools view -Su -b sample_${i}_align.sam | samtools sort - sample_${i}_alns.sorted
rm sample_${i}_align.sam
done

