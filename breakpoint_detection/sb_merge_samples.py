#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:22:52 2020

@author: naokoiida

merge sample files

(1) prepare a list file with sample name list to merge.
example list.txt (id,group) * id might be changed after mapping.
T110,1
T111,1

(2) merge_name and suffix
a role of the file name:
BBNX_{id}_{suffix}.txt 
output file
{merge_name}

"""

def merge_samples(list_file, data_root):    
    import pandas as pd
    
    a_col = ['chr', 'position']
    a_df = pd.DataFrame(index=[], columns=a_col)
    
    with open(list_file, 'r') as hin:
        for raw in hin:
            sample = raw.rstrip('\n').split(',')
            id=sample[0]
            group=sample[1]
            pre=sample[2]
            file_r = f'{data_root}/{group}/{id}/BBNX_{pre}_break.cs.ann.fil.strict.txt'
            
            df = pd.read_csv(file_r, delimiter='\t',usecols=[0,1,4], header=0, dtype={0:'str',1:'str',2:'str'}) 
            res = pd.merge(a_df, df, on=['chr','position'], how="outer")
            a_df = res

    out_file = "./merge/"+group+"_BBNX_break.cs.ann.fil.strict.txt" 
          
    a2_df = pd.concat([a_df,pd.DataFrame(a_df.count(axis=1)-2,columns=['count'])],axis=1)
    a3_df = a2_df.fillna(0)
    
    a3_df.to_csv(out_file, index=False, header=True, sep='\t')

            
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--list", action="store", dest="list",
                            help = "sample number list", required=True) 
    parser.add_argument("--data_root", action="store", dest="data_root",
                            help = "example: ../../data/DDS-study", required=True)

    
    args = parser.parse_args()
    
    list_file = args.list
    data_root = args.data_root
 
    merge_samples(list_file, data_root)

