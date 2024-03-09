#!/bin/bash


#input=./sample/id_list_all.txt
input=$1
python ../../code/breakpoint_detection/sb_proc_list2.py -input_file ${input} 

