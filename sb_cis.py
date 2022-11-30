#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 12:43:09 2019
separate chr position by a sliding window and calculate insertion counts par sample.
in_file = 'merge/organoid_BBNX_break.cs.ann.fil.strict.txt'
[0]chr	
[1]position
[2:-2]sample
[-1]total insertion counts
   
out_file = 'merge/organoid_BBNX_break.cs.ann.fil.strict.cis.txt'

@author: genome
"""

def sb_slide_window(input_file, size, size_name):

    import pathlib
    import pandas as pd
    import os
    
    
    def window_allocate(pos, size):
        import math
        from decimal import Decimal, ROUND_HALF_UP
        
        p1 = math.floor(int(pos)/size)*size
        
        p2a = Decimal(str(int(pos)/size)).quantize(Decimal('0'), rounding=ROUND_HALF_UP)
        p2 = round((float(p2a)-0.5)*size)
        if p2 <0:
            p2 =0
        return p1, p2
#input_file = 'test_strict.txt'    
#size_name='10k'
    in_path = pathlib.Path(input_file)
    name = in_path.stem
    tmp_file = './merge/' + name + '.win' +size_name +'_tmp.txt'
    win_file = './merge/' + name + '.win' +size_name +'.txt'
    input_file_path = './merge/' + input_file
    print(tmp_file)  
    print(input_file_path)    
    with open(input_file_path, 'r') as hin, open(tmp_file, 'w') as tout:
        #columns = ['chr','position','sample1',....,(count)]
        header = hin.readline().rstrip('\n').split('\t')
        new_header = '\t'.join([header[0],'position_start','position_end',]+header[2:])+'\n'
        tout.write(new_header)
        for line in hin:
            F = line.rstrip('\n').split('\t')
            chr = F[0]
            data = '\t'.join(F[2:])
            res_pos = window_allocate(F[1], size)
            end0 = res_pos[0]+size-1 
            end1 = res_pos[1]+size-1 
            row0 = str(chr)+'\t'+str(res_pos[0])+'\t'+str(end0)+'\t'+data+'\n' 
            row1 = str(chr)+'\t'+str(res_pos[1])+'\t'+str(end1)+'\t'+data+'\n'
            tout.write(row0)
            tout.write(row1)
        
        tmp_df = pd.read_csv(tmp_file, header=0, sep='\t')
        group = tmp_df.groupby(['chr','position_start','position_end'], as_index=False)
        agg_group = group.agg(sum) 
        agg_group.to_csv(win_file, index=False, header=True, sep='\t')
    os.remove(tmp_file) 
    
"""    
sb_cal_pvalue        
poisson
from math import exp
from math import factorial
Pr(X=k)=e^λ(λ)^x/x!
prob = (lamb**x)*exp(-lamb)/factorial(x)
"""

def sb_cal_pvalue(input_file, size, size_name):
    
    import pandas as pd
    from scipy.stats import poisson
    import pathlib
    
    in_path = pathlib.Path(input_file)
    name = in_path.stem
    win_file = './merge/' + name + '.win' +size_name +'.txt'
    p_file = './merge/' + name + '.win' +size_name +'.p.txt'
#win_file = './merge/organoid_BBNX_break.cs.ann.fil.strict.win10k.txt'
    df = pd.read_csv(win_file, delimiter='\t', header=0)  
    #total insertion sites
    ins_total = df['count'].sum()
    #genome without N length
    genome_len = 2740108151
    #lamb_b = float(ins_total/(genome_len*sample_n))
    lamb = float((ins_total/genome_len)*size)

    with open(win_file, 'r') as win:
        with open(p_file, 'w') as pout:
            F1 = win.readline().rstrip('\n').split('\t')
            data_index = F1[3:-1]
            rec = F1[0]+'\t'+F1[1]+'\t'+F1[2]+'\t'+F1[-1]+'\tprob,lambda='+str("{:.5E}".format(lamb))+'\tp-value\tsample\n'
            pout.write(rec)
            for line in win:
                F = line.rstrip('\n').split('\t')
                data = F[3:-1]
                n = 0
                data_list =[]
                for i in data:
                    if float(i) >0:
                        data_list.append(data_index[n])
                        data_rec = ','.join(data_list)
                    n +=1
                #ins_count par 10k
                x = float(F[-1])  
                prob = poisson.pmf(x,lamb)
                #cprob = poisson.cdf(x-1,lamb)
                #p = 1 - cprob
                surfun = poisson.sf(x-1, lamb)
                p_float = "{:.2E}".format(surfun)
                rec = F[0]+'\t'+F[1]+'\t'+F[2]+'\t'+F[-1]+'\t'+str(prob)+'\t'+str(p_float)+'\t'+ str(data_rec) +'\n'
                pout.write(rec)
    
def sb_stech_windows(input_file, size_name):
    import pandas as pd
    import pathlib
    import os
#input_file = "organoid2019_BBNX_break.cs.ann.fil.strict.txt"
#size_name = "10k"
    in_path = pathlib.Path(input_file)
    name = in_path.stem
    p_file = './merge/' + name + '.win' +size_name +'.p.txt'
    stech_tmp_file = './merge/' + name + '.win' +size_name +'.p.stech_tmp.txt'
    stech_file = './merge/' + name + '.win' +size_name +'.p.stech.txt'

    df = pd.read_csv(p_file, delimiter='\t', header=0) 
    df_s = df.sort_values(['chr','position_start'])
    df_s.to_csv(stech_tmp_file, index=False, header=True, sep='\t')
    
    num_lines = len(open(stech_tmp_file).readlines())
    n = 0
    with open(stech_tmp_file, 'r') as win, open(stech_file, 'w') as sout:
        #read the header
        header = win.readline().rstrip('\n').split('\t')
        n +=1
        new_header = '\t'.join([header[0],header[1],header[2],header[3],header[5]]) +'\t'+'min_p-value'+'\t'+header[6]+'\n'
        sout.write(new_header)
        #read the first line
        F1 = win.readline().rstrip('\n').split('\t')
        n +=1
        chrom1 = F1[0]
        start0 = F1[1]
        start1 = F1[1]
        end1 = F1[2]
        count1 = F1[3]
        p1 = float(F1[5])
        
        count_list =[count1]
        p_list =[p1]
        
        sample_list = [i for i in F1[6].split(',')]
        
        for line in win:
            n +=1
            F = line.rstrip('\n').split('\t')
            chrom = F[0]
            start = F[1]
            end = F[2]
            count = F[3]
            p = float(F[5])
            sample = [i for i in F[6].split(',')]
            
            if chrom == chrom1 and int(start1)<int(start) and int(start)<int(end1):
                
                count_list.append(count)
                p_list.append(p)
                sample_list = sample_list + sample
                chrom1 = chrom
                start1 = start
                end1 = end
            
            else:
                
                count_set = ','.join(count_list)
                p_set_s = [str(s) for s in p_list]
                p_set = ','.join(p_set_s)
                p_min = min(p_list)
                sample_set = ','.join(set(sample_list))
                rec = str(chrom)+'\t'+str(start0)+'\t'+str(end1)+'\t'+str(count_set)+'\t'+str(p_set)+'\t'+str(p_min)+'\t'+str(sample_set)+'\n'
                sout.write(rec)
                chrom1 = chrom
                start0 = start
                start1 = start
                end1 = end
                
                count_list =[count]
                p_list =[p]
                sample_list = sample
                
            if n == num_lines:
                count_set = ','.join(count_list)
                p_set_s = [str(s) for s in p_list]
                p_set = ','.join(p_set_s)
                p_min = min(p_list)
                sample_set = ','.join(set(sample_list))
                rec = str(chrom)+'\t'+str(start0)+'\t'+str(end1)+'\t'+str(count_set)+'\t'+str(p_set)+'\t'+str(p_min)+'\t'+str(sample_set)+'\n'
                sout.write(rec)
            
    os.remove(stech_tmp_file) 

"""

"""


if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--infile",action="store", metavar = "folder name", default = None, type = str,
                            help = "folder name with input files. BBNXcombine_*_onebreak.bp.classify.gene.fil.txt") 
    parser.add_argument("--size", action="store", metavar = "window size", default = 10000, type = int,
                            help = "window size") 
    parser.add_argument("--size_name",action="store", metavar = "tag for file name", default = '10k', type = str,
                            help = "tag for file name. [10k]") 
            
    args = parser.parse_args()
    
    input_file = args.infile
    size = args.size
    size_name = args.size_name
 
    sb_slide_window(input_file, size, size_name)
    sb_cal_pvalue(input_file, size, size_name)
    sb_stech_windows(input_file, size_name)
    
    #sb_stech_windows("test.txt","10k")
 
