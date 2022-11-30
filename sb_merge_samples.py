#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 15:22:52 2020

@author: naokoiida

merge sample files

(1) prepare a list file with sample name list to merge.
example list.txt
110
111

(2) merge_name and suffix
a role of the file name:
./break/BBNX_{sample[0]}_{suffix}.txt 
output file
./merge/{merge_name}_{suffix}.txt'

"""

def merge_samples(list_file, suffix, merge_name):    
    import pandas as pd
    
    out_file = f'./merge/{merge_name}_BBNX_{suffix}.txt'
    
    a_col = ['chr', 'position']
    a_df = pd.DataFrame(index=[], columns=a_col)
    
    with open(list_file, 'r') as hin:
        for raw in hin:
            sample = raw.rstrip('\n').split('\t')
            file_r = f'./break/BBNX_{sample[0]}_{suffix}.txt'
            
            df = pd.read_csv(file_r, delimiter='\t',usecols=[0,1,4], header=0, dtype={0:'str',1:'str',2:'str'}) 
            res = pd.merge(a_df, df, on=['chr','position'], how="outer")
            a_df = res
            
        a2_df = pd.concat([a_df,pd.DataFrame(a_df.count(axis=1)-2,columns=['count'])],axis=1)
        a3_df = a2_df.fillna(0)
    
        a3_df.to_csv(out_file, index=False, header=True, sep='\t')

            
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--list", action="store", dest="list",
                            help = "sample number list", required=True) 
    parser.add_argument("--suffix", action="store", dest="suffix",
                            help = "BBNX_{sample[0]}_{suffix}.txt", required=True)
    parser.add_argument("--merge_name", action="store", dest="merge_name",
                            help = "merge name such as tissue and tumor", required=True) 
    

    args = parser.parse_args()
    
    list_file = args.list
    suffix = args.suffix
    merge_name = args.merge_name
 
    merge_samples(list_file, suffix, merge_name)

