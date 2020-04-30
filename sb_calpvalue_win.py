#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 09:50:21 2019

@author: genome
poisson
from math import exp
from math import factorial
Pr(X=k)=e^λ(λ)^x/x!
prob = (lamb**x)*exp(-lamb)/factorial(x)
SB_calvalue.py
"""
def sb_cal_pvalue(input_dir, size, name):
    
    import pandas as pd
    from scipy.stats import poisson
    
    in_file = f'./{input_dir}/window_{name}/{input_dir}_BBNXcombine_all_onebreak.count_{name}.txt'
    out_file = f'./{input_dir}/window_{name}/{input_dir}_BBNXcombine_all_onebreak.count_{name}.p.txt'

    df = pd.read_csv(in_file, delimiter='\t', header=0)
    col = df.columns.tolist()
    sample_col = [l for l in col if l.startswith('BBNX')]    
    #sample_n = len(sample_col)
    ins_total = df['sum'].sum()
    #genome without N length
    genome_len = 2740108151
    #lamb_b = float(ins_total/(genome_len*sample_n))
    lamb = float((ins_total/genome_len)*size)


    with open(in_file, 'r') as in1:
        with open(out_file, 'w') as out1:
            for line in in1:
                F = line.rstrip('\n').split('\t')
                l = line.rstrip('\n')
                if F[1] == 'start' and F[-1] == 'sum': #header check
                    rec = 'chr\tstart\tend\tsum\tlambda=(ins_total/genome_len)*window_size\tprob\tp-value\t' + l +'\n'
                    out1.write(rec)
                elif F[-1] == 'NaN':
                    rec = F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[-1] + '\tNaN\tNaN\tNaN\t' + l + '\n'
                    out1.write(rec)
                else:
                    x = float(F[-1])  #ins_count par gene
                    prob = poisson.pmf(x,lamb)
                    #cprob = poisson.cdf(x-1,lamb)
                    #p = 1 - cprob
                    surfun = poisson.sf(x-1, lamb)
                    p_float = "{:.2E}".format(surfun)
                    rec = F[0] + '\t' + F[1] + '\t' + F[2] + '\t' + F[-1] +'\t' + str(lamb) + '\t' + str(prob) + '\t' + str(p_float) + '\t' + l + '\n'
                    out1.write(rec)

    #add gene info

    gene_file = f'./{input_dir}/window_{name}/{input_dir}_BBNXcombine_all_onebreak.bp.classify.gene.fil.geneinfo.{name}.txt'                        
    in_file = f'./{input_dir}/window_{name}/{input_dir}_BBNXcombine_all_onebreak.count_{name}.p.txt'
    out_file = f'./{input_dir}/window_{name}/{input_dir}_BBNXcombine_all_onebreak.count_{name}.p.gene.txt'
    
    df_gene = pd.read_csv(gene_file, header=0, sep='\t', dtype=str)
    #df_gene.columns = ['chr','start','end','break_info']
    df_gene2 = df_gene.set_index(["chr","start","end"], drop=True)
    
    df = pd.read_csv(in_file, header=0, sep='\t', dtype=str) #,index=['chr','start','end','position','gene','insertion_count']
    df2 = df.set_index(["chr","start","end"], drop=True)
    
    res = pd.merge(df_gene2, df2, left_index=True, right_index=True, how="right")
    res2 = res.fillna(0)
    xrec = res2.reset_index()
    xrec.to_csv(out_file, index=False, header=True, sep='\t')


if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("input_dir", metavar = "folder name", default = None, type = str,
                            help = "folder name with input files. BBNXcombine_*_onebreak.bp.classify.gene.fil.txt") 
    parser.add_argument("size", metavar = "window size", default = None, type = str,
                            help = "window size") 
    parser.add_argument("name", metavar = "tag for file name", default = None, type = str,
                            help = "tag for file name. [10k]") 
            
    args = parser.parse_args()
    
    input_dir = args.input_dir
    size = args.size
    name = args.name
 
    sb_cal_pvalue(input_dir, size, name)

#https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.poisson.html
sb_cal_pvalue('Primary', 10000, '10k')