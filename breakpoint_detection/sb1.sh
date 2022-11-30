#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=8G,s_vmem=8G
#$ -pe def_slot 2

ID=$1
ENZ=$2
GROUP=$3

DIR_PATH=../../data/DDS-study
BAM=${DIR_PATH}/${GROUP}/${ID}/${ENZ}_T${ID}.bam
PREFIX=${DIR_PATH}/${GROUP}/${ID}/${ENZ}_T${ID}

python ../../code/breakpoint_detection/breakpoint_detectorS.py --infile ${BAM} --outfile ${PREFIX}_break.txt --re ${ENZ}

python ../../code/breakpoint_detection/consensus_maker.py --infile ${PREFIX}_break.txt --outfile ${PREFIX}_break.cs.txt --cut_off 0.8

python ../../code/breakpoint_detection/annot_sbseq.py --input_file ${PREFIX}_break.cs.txt --output_file ${PREFIX}_break.cs.ann.fil.strict.txt

