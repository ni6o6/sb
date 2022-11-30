#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:43:17 2019

@author: genome

after filtering.
"""

'''
file list:
    NX005.fastq.gz	<\t> NX_005<\t> NX or BB
    
BB_{i}
NX_{i}
input_file header
[0]chr	
[1]position(0-start)	
[2]ori_position    
[3]family_size  
[4]sb_direction	
[5]adj_seq
[6]soft_clip_len   
[7]sw_match        
[8]sw_match_ratio  
[9]sb-seq| 
[10]|genome2bases	


output_file
[0] chr
[1] position
[2] sb_direction
[3] sb_genome_seq
[4] BBNX_ID
'''

def sb_combBBNX(id, dir_path):
    import pandas as pd
    import os
    
    bb_file = f'{dir_path}/BB_{id}_break.cs.ann.fil.strict.txt'
    nx_file = f'{dir_path}/NX_{id}_break.cs.ann.fil.strict.txt'
    
    out_tmp = f'{dir_path}/BBNX_{id}_break.cs.ann.fil.strict_tmp.txt'
    out_file = f'{dir_path}/BBNX_{id}_break.cs.ann.fil.strict.txt'

    if os.path.isfile(bb_file):
        bb_data = pd.read_csv(bb_file, delimiter='\t',usecols=[0,1,4,5], header=0, dtype={0:'str',1:'int',2:'str',3:'str'}) 
        bb_data.columns = ['chr', 'position','bb_sb_direction', 'bb_seq']
        
    else:
        bb_col = ['chr', 'position','bb_sb_direction', 'bb_seq']
        bb_data = pd.DataFrame({'chr':"c0", 'position':1,'bb_sb_direction':"c2",'bb_seq':"c3"},index=[1], columns=bb_col)
        
    if os.path.isfile(nx_file):
        nx_data = pd.read_csv(nx_file, delimiter='\t',usecols=[0,1,4,5], header=0, dtype={0:'str',1:'int',2:'str',3:'str'}) 
        nx_data.columns = ['chr', 'position','nx_sb_direction', 'nx_seq']

    else:
        nx_col = ['chr', 'position','nx_sb_direction', 'nx_seq']
        nx_data = pd.DataFrame({'chr':"c0", 'position':1,'bb_sb_direction':"c2",'bb_seq':"c3"},index=[1], columns=nx_col)

    bbnx_data = pd.merge(bb_data, nx_data, on=['chr', 'position'], how='outer', indicator=True)
    bbnx_data.to_csv(out_tmp, index=False, header=True, sep='\t')
    
    with open(out_tmp, 'r') as hin, open(out_file, 'w') as hout:
        next(hin)
        col = f'BBNX_{id}'
        header = '\t'.join(['chr','position','sb_direction','sb_genome_seq',col])
        hout.write(header+'\n')
        
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] == "c0": continue
            
            if F[6] == 'right_only':
                rec = str(F[0])+'\t'+str(F[1])+'\t'+str(F[4])+'\t'+str(F[5])+'\t1\n' 
            else:
                rec = str(F[0])+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(F[3])+'\t1\n'
            hout.write(rec)
            
    os.remove(out_tmp)
       
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    parser.add_argument("--id", action="store", dest="id",
                            help = "sample name that is common for BB and NX.", required=True)
    parser.add_argument("--dir_path", action="store", dest="dir_path",
                            help = "./break to data", required=True)

    args = parser.parse_args()
    
    id = args.id
    dir_path = args.dir_path
 
    sb_combBBNX(id, dir_path)
    


'''

sample='103'
suffix='_break.cs.ann.fil.strict'
dir_path='./break/'
'''



        
