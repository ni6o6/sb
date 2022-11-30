#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=4G,s_vmem=4G
#$ -pe def_slot 2

data_dir=../../data
file=fastq_list.txt

for line in $(cat ${file}); do
    input_id=`echo ${line} | cut -d',' -f1`
    echo ${input_id}
    output_id=`echo ${line} | cut -d',' -f2`
    enz=`echo ${line} | cut -d',' -f3`
    group=`echo ${line} | cut -d',' -f4`

    input=../../fastq/s${input_id}.fastq.gz
    output_prefix=${data_dir}/DDS-study/${group}/${output_id}/${enz}_T${output_id}

    mkdir -p ${data_dir}/DDS-study/${group}/${output_id}/
    bash ../../code/bwa/smap.sh ${input} ${output_prefix}

done

