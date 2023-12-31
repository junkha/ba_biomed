#! /usr/bin/bash

snvs_with_segments=/home/junkhann/daten/snvs_with_segments/$1
onek1k=/home/junkhann/daten/onek1k/onek1k_subset_columns.tsv
dest_extra_columns=/home/junkhann/daten/snvs_with_segments_extra_columns
dest_joined=/home/junkhann/daten/joined_mmml_onek1k

for file in $snvs_with_segments/*; do
    pid=$(echo $file | tr -dc '0-9')

    echo "Add columns on PID ${pid}"
    python3 add_columns.py -f $file -o ${dest_extra_columns}/extra_columns_${pid}.vcf.gz
    
    echo "Intersect on PID ${pid}"
    intersectBed -wa -wb -a ${dest_extra_columns}/extra_columns_${pid}.vcf.gz -b $onek1k | gzip > ${dest_joined}/mmmml_onek1k_${pid}.vcf.gz

    echo "Done for PID ${pid}"
done