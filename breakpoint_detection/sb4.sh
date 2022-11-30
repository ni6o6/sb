#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=8G,s_vmem=8G
#$ -pe def_slot 2

id_list=$1
group_name=$2
code_path=/home/naiida/SleepingBeauty/code/breakpoint_detection
data_root=../../data/DDS-study
simple=../../reference/mm10_simpleRepeat_modsort.bed.gz 
segment=../../reference/mm10_segmentalDups_modsort.bed.gz

python ${code_path}/sb_merge_samples.py --list ${id_list} --data_root ${data_root} 

python ${code_path}/sb_cis.py --infile ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.txt

python ${code_path}/add_mm10repeats.py ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.stech.txt ${simple} ${segment}

input_file=./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.stech.Simprep.Segdups.txt
path_to_db=../../reference/
python ${code_path}/sb_postproc_backTo1.py --input_file ${input_file} --path_to_data_dir ${data_root} --path_to_db ${path_to_db} --group ${group_name} 

rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.stech.Simprep.Segdups.ins1base.txt
rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.stech.Simprep.Segdups.txt
rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.stech.txt
rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.p.txt
rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.win10k.txt
rm ./merge/${group_name}_BBNX_break.cs.ann.fil.strict.txt
