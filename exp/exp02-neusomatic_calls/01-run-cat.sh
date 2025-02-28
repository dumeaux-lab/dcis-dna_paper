#!/bin/bash
filename="$1"

while read -r field1 field2 field3 field4 ;  
do
echo $field1
cat <(cat /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/*/Wrapper/Ensemble.s*.tsv |grep CHROM|head -1) \
    <(cat /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/*/Wrapper/Ensemble.s*.tsv |grep -v CHROM) | sed "s/nan/0/g" > /path/to/bc_dcis/DNAseq_allBatches/neusomatic2/$field1/ensemble_ann.tsv
    
done < "$filename"

