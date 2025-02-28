#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1

python /path/to/tools/neusomatic/neusomatic/python/call.py \
--candidates_tsv /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_precall_ensemble/dataset/*/candidates*.tsv \
--reference /path/to/db/gatk_hg38/Homo_sapiens_assembly38.fasta \
--out /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/work_call_ensemble \
--checkpoint /path/to/tools/neusomatic/neusomatic/models/NeuSomatic_v0.1.4_ensemble_SEQC-WGS-Spike.pth \
--num_threads 10 \
--ensemble \
--batch_size 100

done < "$filename"

