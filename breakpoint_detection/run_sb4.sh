#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=4G,s_vmem=4G
#$ -pe def_slot 1


GLIST="1 11 12 13 14 15 16 2 21 22 23 24 2D 2P 2PM 3 31 32 33 34 35 4 5 6 Small"
#GLIST="31 Small"
#GLIST="6"

for GROUP in ${GLIST}; do

code_path=/home/naiida/SleepingBeauty/code/breakpoint_detection
qsub ${code_path}/sb4.sh ./sample/id_list_${GROUP}.txt ${GROUP}

done

