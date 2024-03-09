#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Naoko Iida

file list:
[0] id 
[1] Group
[2] subgroup(if any or same with Group)
"""

def sb_proc_list(input_file):

        
    with open(input_file, 'r') as hin:
        id_dict = {}
        
        for line in hin:
            F = line.rstrip('\n').split(',')
            id = F[0]
            if id in id_dict.values(): continue
            else:
                # id,group,id for file prefix
                if id in ["Villin149ST1","Villin149ST2","Villin149ST3", "Villin202PT1", "Villin202PT2"]:
                    id_dict.setdefault(F[2],[]).append(F[0] +","+ F[1] +","+ F[0])
                else:
                    id_dict.setdefault(F[2],[]).append(F[0] +","+ F[1] +","+ "T" + F[0])                    

    for key in id_dict:
        output_file = "./sample/id_list_"+str(key)+".txt"        
        with open(output_file, 'w') as hout:           
            value_list=id_dict.get(key)
            rec = '\n'.join(value_list)
            hout.write(rec+"\n")

       
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    parser.add_argument("-input_file", action="store", dest="input_file",
                            help = "input file", required=True)

    args = parser.parse_args()
     
    sb_proc_list(args.input_file)
    





        