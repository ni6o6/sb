#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Naoko Iida

file list:
[0] seq id
[1] id 
[2] NX/BB 
[3] Group
"""

def sb_proc_list(input_file, output_file):

        
    with open(input_file, 'r') as hin:
        id_dict = {}
        id_list = []
        for line in hin:
            F = line.rstrip('\n').split(',')
            if F[1] in id_list: continue
            else:
                id_dict[F[1]] = F[3]
            
    with open(output_file, 'w') as hout:           
        for id, group in id_dict.items():
            hout.write(str(id)+","+str(group)+"\n")

       
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    parser.add_argument("-input_file", action="store", dest="input_file",
                            help = "input file", required=True)
    parser.add_argument("-output_file", action="store", dest="output_file",
                            help = "unified id list", required=True)

    args = parser.parse_args()
     
    sb_proc_list(args.input_file, args.output_file)
    





        