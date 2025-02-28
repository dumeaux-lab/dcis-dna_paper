#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1
echo $field2
echo $field4

# EJM note: Only changed path to files, all other code left unchanged
python /path/to/tools/neusomatic/neusomatic/python/preprocess.py \
	--mode call \
	--reference /path/to/db/gatk_hg38/Homo_sapiens_assembly38.fasta \
	--region_bed /path/to/db/GRCh38/SureSelectHumanExonsv7/S31285117_hs_hg38/S31285117_Regions.bed \
	--tumor_bam $field2 \
	--normal_bam $field4 \
	--work /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_precall_ensemble \
	--ensemble_tsv /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/ensemble_ann.tsv \
	--min_mapq 10 \
	--num_threads 10 \
	--scan_alignments_binary /path/to/tools/neusomatic/neusomatic/bin/scan_alignments


done < "$filename"

