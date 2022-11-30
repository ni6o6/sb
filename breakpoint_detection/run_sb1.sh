#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=4G,s_vmem=4G
#$ -pe def_slot 1

file=$1

for line in $(cat ${file}); do
    id=`echo ${line} | cut -d',' -f2`
    echo ${id}
    enz=`echo ${line} | cut -d',' -f3`
    group=`echo ${line} | cut -d',' -f4`
    
    bash ../../code/breakpoint_detection/sb1.sh ${id} ${enz} ${group}

done
