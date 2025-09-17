#!/bin/bash

input_file=$(realpath $1)
output_file=$(realpath $2)


for file in ${input_file}/*_pph.tsv; do

    echo $file
    
    if [ -f $output_file ]; then
        awk -F "\t" 'NR>1 && ($6>=0.7 || $7>=0.7)' $file >>$output_file
    else
        echo -e "eqtl_data_from\tgene_id\tPP.H0\tPP.H1\tPP.H2\tPP.H3\tPP.H4" >$output_file
        awk -F "\t" 'NR>1 && ($6>=0.7 || $7>=0.7)' $file >>$output_file
    
fi
done






