# rnaseq_aws_benchmarking
A roadmap for setting up and benchmarking RNA-seq alignment tools on AWS.

Tutorial: [RNA-seq Alignment Benchmarking Through AWS](https://github.com/Slandert/rnaseq_aws_benchmarking/wiki)

Code workflow:

Setup and alignment of benchmark dataset reads:

UnixScripts:
1. Get_benchmark_fastas.sh
2. STAR_align_reads.sh
3. rsem_calculate_expression.sh
4. consolidate_rsem.sh
5. kallisto_quant.sh
6. salmon_quant.sh

Comparison of aligner output:

compareAlignments/compareAlignments.R

