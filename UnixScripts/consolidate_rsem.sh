for i in {1..30};
do
mkdir download-data/sample${i}
cp sample_${i}/Quant.genes.results download-data/sample${i}
cp sample_${i}/Quant.isoforms.results download-data/sample${i}

done
