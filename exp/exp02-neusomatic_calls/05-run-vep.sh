#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1

/path/to/tools/ensembl-vep/vep --cache --offline --format vcf --vcf --force_overwrite --fork 10 \
-i /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble/NeuSomatic_ensemble.vcf \
-o /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble/NeuSomatic_ensemble_annot.vcf

done < "$filename"
