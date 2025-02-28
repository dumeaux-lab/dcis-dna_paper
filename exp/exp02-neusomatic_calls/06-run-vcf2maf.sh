#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1

perl /path/to/tools/mskcc-vcf2maf-5453f80/vcf2maf.pl \
--vep-path /path/to/tools/ensembl-vep \
--vep-data /path/to/.vep \
--input-vcf /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble/NeuSomatic_ensemble_annot.vcf \
--output-maf /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_postcall_ensemble/NeuSomatic_ensemble_annot.maf \
--tumor-id $field1 \
--vcf-tumor-id SAMPLE \
--species homo_sapiens --ncbi-build GRCh38 \
--retain-info SCORE --retain-fmt AF \
--ref-fasta /path/to/db/gatk_hg38/Homo_sapiens_assembly38.fasta

done < "$filename"
