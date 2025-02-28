#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1
echo $field2

python /path/to/tools/neusomatic/neusomatic/python/postprocess.py \
--reference /path/to/db/gatk_hg38/Homo_sapiens_assembly38.fasta \
--tumor_bam $field2 \
--pred_vcf /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_call_ensemble/pred.vcf \
--candidates_vcf /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_precall_ensemble/work_tumor/filtered_candidates.vcf \
--ensemble_tsv /path/to/bc_dcis/DDNAseq_allBatches/neusomatic2/$field1/ensemble_ann.tsv \
--output_vcf /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble/NeuSomatic_ensemble.vcf \
--work /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble


done < "$filename"

