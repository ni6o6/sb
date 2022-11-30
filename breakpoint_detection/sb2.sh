#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=4G,s_vmem=4G
#$ -pe def_slot 2

ID=$1
GROUP=$2

ROOT_PATH=../../data/DDS-study
DIR_PATH=${ROOT_PATH}/${GROUP}/${ID}

python ../../code/breakpoint_detection/sb_combBBNX_break.py --id T${ID} --dir_path ${DIR_PATH}

python ../../code/breakpoint_detection/sb_add_gene.py --input_file ${DIR_PATH}/BBNX_T${ID}_break.cs.ann.fil.strict.txt --output_file ${DIR_PATH}/BBNX_T${ID}_break.cs.ann.fil.strict.addgene.txt --gene_file ../../reference/gene.bed.gz

