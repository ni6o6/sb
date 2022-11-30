#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=8G,s_vmem=8G
#$ -pe def_slot 2

input=./sample/mapping_list.txt
output=./sample/mapping_list_uniq.txt

python ../../code/breakpoint_detection/sb_proc_list.py -input_file ${input} -output_file ${output}

for line in $(cat ${output}); do
    id=`echo ${line} | cut -d',' -f1`
    echo ${id}
    group=`echo ${line} | cut -d',' -f2`
    
    bash ../../code/breakpoint_detection/sb2.sh ${id} ${group}

done
