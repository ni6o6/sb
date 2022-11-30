#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:31:48 2020

@author: genome
python mm10_add_repeats.py input_file, simple_file, segdup_file
"""

def Repeat_filter(input_file, simple_file, segdup_file):
    
    import pysam
    import pathlib
    import os
    
    simple_tb = pysam.TabixFile(simple_file)
    segdup_tb = pysam.TabixFile(segdup_file)
    
    in_path = pathlib.Path(input_file)
    out_file1 = "./merge/"+in_path.stem +'.Simprep.txt'
    out_file2 = "./merge/"+in_path.stem +'.Simprep.Segdups.txt'
  
    h1out = open(out_file1, 'w') 
    
    target = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
    
    with open(input_file, 'r') as h1in:
        #next(hin)
        header1 = h1in.readline().rstrip('\n')
        out1_header = header1 + '\tsimpleRepeat_score\tsimpleRepeat_seq'
        print(out1_header, file = h1out)
    
        for line in h1in:
            ln = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
            
            tchr = F[0]
            tpos_s = F[1]
            tpos_e = F[2]
    
            if tchr in target:
                records = simple_tb.fetch(tchr, int(tpos_s)-1, int(tpos_e))
                
                simple_score = '0'
                simple_seq = '-'
                if records is not None:
                    for records_line in records:
                        R = str(records_line).split('\t') 
                        simple_score = str(R[3])
                        simple_seq = str(R[4])
                result = ln +'\t'+ simple_score +'\t'+simple_seq + '\n'
                h1out.write(result)
            
    h1out.close()
    
    
    #filter out segmental dups
    with open(out_file1, 'r') as h2in, open(out_file2, 'w') as h2out:
        #next(hin)
        header2 = h2in.readline().rstrip('\n')
        out2_header = header2 + '\tsegmentalDups'
        print(out2_header, file = h2out)
    
        for line in h2in:
            ln = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
    
            tchr = F[0]
            tpos_s = F[1]
            tpos_e = F[2]
           
            if tchr in target:
                
                records = segdup_tb.fetch(tchr, int(tpos_s)-1, int(tpos_e))
                
                segdup_rec = '-'
                if records is not None: 
                    for records_line in records:
                        R = str(records_line).split('\t') 
                        segdup_rec = str(R[3])
                
                #if segdup_rec == '-':
                result = ln +'\t'+ segdup_rec + '\n'
                h2out.write(result)
                
    h2out.close()
    
    os.remove(out_file1)

    

           
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input_file", metavar = "input_file", default = None, type = str,
                            help = "Path to input file") 
        
    parser.add_argument("simple_file", metavar = "simple_file", default = "hg38_simpleRepeat.bed.gz", type = str,
                            help = "Path to simple repeats data file") 
    
    parser.add_argument("segdup_file", metavar = "segdup_file", default = "hg38_simpleRepeat.bed.gz", type = str,
                            help = "Path to simple repeats data file")     
        
    args = parser.parse_args()
    
    input_file = args.input_file
    simple_file = args.simple_file
    segdup_file = args.segdup_file
    
    Repeat_filter(input_file, simple_file, segdup_file)

"""         
input_file='primary1d_BBNX_break.cs.ann.fil.strict.win10k.p.stech.addgene.txt'
simple_file='/Volumes/NIIDA_SleepingBeauty/database/mm10_simpleRepeat_modsort.bed.gz'
segdup_file='/Volumes/NIIDA_SleepingBeauty/database/mm10_segmentalDups_modsort.bed.gz'
Repeat_filter(input_file, simple_file, segdup_file)
"""
