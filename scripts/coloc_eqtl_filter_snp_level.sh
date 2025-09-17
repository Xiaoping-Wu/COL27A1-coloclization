#!/bin/bash


input_file=$(realpath $1)
output_file=$(realpath $2)


for file in ${input_file}/*_snp_level.tsv.gz; do

    echo $file
    
    if [ -f $output_file ]; then
        zcat $file |awk -F  "\t" 'BEGIN{OFS="\t"} NR>1 && $11>0.7 {print}'  >>$output_file
        
    else
        echo -e "snp\tV.df1\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tgene_id\teqtl_data_from" >$output_file
        zcat $file |awk -F  "\t" 'BEGIN{OFS="\t"} NR>1 &&  $11>0.7 {print}'  >>$output_file
fi
done





#input_file=$(realpath $1)
#snplist=$(realpath $2)
#output_file=$(realpath $3)


#for file in ${input_file}/*_snp_level.tsv.gz; do

#    echo $file
    
#    if [ -f $output_file ]; then
#        zgrep -F -w -f $snplist $file |awk -F  "\t" 'BEGIN{OFS="\t"} $11>0.7 {print}'  >>$output_file
        
#    else
#        echo -e "snp\tV.df1\tz.df1\tr.df1\tlABF.df1\tV.df2\tz.df2\tr.df2\tlABF.df2\tinternal.sum.lABF\tSNP.PP.H4\tgene_id\teqtl_data_from" >$output_file
#        zgrep -F -w -f $snplist $file |awk -F  "\t" 'BEGIN{OFS="\t"} $11>0.7 {print}'  >>$output_file
#fi
#done






