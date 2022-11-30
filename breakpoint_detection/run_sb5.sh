#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=8G,s_vmem=8G
#$ -pe def_slot 2

GLIST="1 11 12 13 14 15 16 2 21 22 23 24 2D 2P 2PM 3 31 32 33 34 35 4 5 6 Small"
#GLIST="31 Small"
#GLIST="6"

for GROUP in ${GLIST}; do

code_path=/home/naiida/SleepingBeauty/code/breakpoint_detection
data=/home/naiida/SleepingBeauty/data/DDS-study
python ${code_path}/sb_makebedfile.py --list ./sample/id_list_${GROUP}.txt --output_file ./bed/${GROUP}_sb_position.bed --path_to_data_dir ${data}

done
