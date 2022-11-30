#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:43:09 2019
@author: naoko iida

"""

def sb_makebedfile(list, output_file, path_to_data_dir):
    
    import os
    
    tmp_file = f'{output_file}_tmp.bed'
    
    with open(list, 'r') as hin, open(tmp_file, 'w') as tmp:
        for raw in hin:
            R = raw.rstrip('\n').split(',')
            id=R[0]
            pre=R[2]
            file_r = path_to_data_dir+"/"+R[1]+"/"+id+"/BBNX_"+pre+"_break.cs.ann.fil.strict.txt"
            with open(file_r, 'r') as rin:
                
                line = rin.readline().rstrip('\n').split('\t')
                name = line[4]
                
                for line in rin:
                    F = line.rstrip('\n').split('\t')
                    pos_end = int(F[1])+1
                    rec = str(F[0])+'\t'+str(F[1])+'\t'+str(pos_end)+'\t'+str(name)+'\t'+str(F[2])+'\n'
                    tmp.write(rec)


    with open(tmp_file, 'r') as jin, open(output_file, 'w') as out:
        for line in jin:
    
            F = line.rstrip('\n').split('\t')
            if F[4] == '+':
                
                rec = str(F[0])+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(F[3])+'\t'+str(1)+'\n'
                out.write(rec)
            else:
                
                rec = str(F[0])+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(F[3])+'\t'+str(-1)+'\n'
                out.write(rec)
    
    os.remove(tmp_file)
            
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--list",action="store", metavar = "list of sample name", default = None, type = str,
                            help = "sample number list file.") 
 
    parser.add_argument("--output_file",action="store", metavar = "output_file", default = 'sample.bed', type = str,
                            help = "output file name") 
    parser.add_argument("--path_to_data_dir", action="store", metavar = "path_to_data_dir", default = None, type = str,
                            help = '/path/to/folder/ for /group/id/BBNX_id_break.cs.ann.fil.strict.txt') 
    args = parser.parse_args()
    
    list = args.list
    output_file = args.output_file
    path_to_data_dir = args.path_to_data_dir
 
    sb_makebedfile(list, output_file, path_to_data_dir)
    



    
