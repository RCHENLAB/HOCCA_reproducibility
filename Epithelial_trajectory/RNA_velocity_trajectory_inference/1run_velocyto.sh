#!/bin/sh
source ~/.condainit
conda activate veloc 
module load samtools/1.10
repeat=/dfs3b/ruic20_lab/junw42/reference/ucsc_gtf/hg38_rmsk.gtf
gtf=/dfs3b/ruic20_lab/junw42/reference/cellranger/refdata-gex-GRCh38-2020-A/genes/genes.gtf

out=$1

velocyto run -b ${out}/filtered_feature_bc_matrix/barcodes.tsv.gz -@ 4 -o $out -m $repeat ${out}/possorted_genome_bam.bam $gtf



