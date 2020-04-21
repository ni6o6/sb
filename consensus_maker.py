#!/usr/bin/env python
"""
Breakpoint Detection
By Naoko Iida
2020.03.26


Inputs: 
	A file with break position and sequence.

Outputs:
	1. a sorted file of input file.
    2. consensus seq of reads with the same chr and break position.
    containing a raw seq around the break position, watson strand seq, sb direction.

usage: breakpoint_detector.py [-h] [--infile INFILE] 
						[--outfile OUTFILE] [--cut_off]
						
optional arguments:
	-h, --help            show this help message and exit
	--infile INFILE       input file
	--outfile OUTFILE     output file
    --cut_off N           default 0.8
"""

from argparse import ArgumentParser
import pandas as pd
import os

def consensus_maker(grouped_reads_list, cut_off):
	# The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no 
	# nucleotide majority reaches above the minimum theshold (--cut_off), the position is considered undefined and an 'N' 
	# is placed at that position in the read.
    read_length = min([len(a) for a in grouped_reads_list])
    nuc_identity_list = [0, 0, 0, 0, 0, 0]  # In the order of T, C, G, A, N, Total
    nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
    consensus_read = ''
    
    # read j1: i1 i2 aaaaaaaaa
    # read j2: aaaataaaaa
    for i in range(read_length):  # Count the types of nucleotides at a position in a read. i is the nucleotide index 
		# within a read in grouped_reads_list
        for j in range(len(grouped_reads_list)):  # Do this for every read that comprises a SMI group. j is the read 
			# index within grouped_reads_list
            try:
                if grouped_reads_list[j][i] == 'T':
                    nuc_identity_list[0] += 1
                elif grouped_reads_list[j][i] == 'C':
                    nuc_identity_list[1] += 1
                elif grouped_reads_list[j][i] == 'G':
                    nuc_identity_list[2] += 1
                elif grouped_reads_list[j][i] == 'A':
                    nuc_identity_list[3] += 1
                elif grouped_reads_list[j][i] == 'N':
                    nuc_identity_list[4] += 1
                else:
                    nuc_identity_list[4] += 1
                nuc_identity_list[5] += 1
            except:
                break
        try:
            for j in [0, 1, 2, 3, 4]:
                if float(nuc_identity_list[j])/float(nuc_identity_list[5]) > cut_off:
                    consensus_read += nuc_key_dict[j]
                    break
                elif j == 4:
                    consensus_read += 'N'
        except:
            consensus_read += 'N'
        nuc_identity_list = [0, 0, 0, 0, 0, 0]  # Reset for the next nucleotide position
    return consensus_read, len(grouped_reads_list)


def run_consensus_maker(infile, outfile, cut_off):
    #reverse complement
    trans = str.maketrans('ATGCatgcN', 'TACGTAGCN')
    
    #sort data
    df = pd.read_csv(f'./break/{infile}', sep ='\t', header=None)
    df.columns = ['chr', 'position(0-start)', 'sb_length', 'read_direction', 'raw_seq', 'watson_seq5to3', 'sb_direction']
    df['i_col'] = df['chr'].str.split('chr', expand=True)[1]
    df_sort = df.sort_values(['i_col','position(0-start)'], ascending=True)
    df_rec = df_sort.drop('i_col',axis=1)
    tmp_file = f'./break/{infile}_sort.txt'
    df_rec.to_csv(tmp_file, sep ='\t', index=None)
    
    grouped_reads_list = []
    position_0 = ''
    
    num_lines = len(open(tmp_file).readlines())
    n = 0
    
    with open(tmp_file, 'r') as hin, open(f'./break/{outfile}', 'w') as hout:
        header = '\t'.join(['chr', 'position(0-start)', 'sb_length', 'read_direction', 'reads', 'raw_seq', 'watson_seq5to3', 'sb_direction'])
        hout.write(str(header)+'\n')
        next(hin)
        #read the first read
        F1 = hin.readline().rstrip('\n').split('\t')
        n +=1
        chr_0 = F1[0]
        position_0 = F1[1]
        chr_position_0 = F1[0]+':'+F1[1]
        sb_length_0 = F1[2]
        read_direction_0 = F1[3]
        sb_direction_0 = F1[6]
        grouped_reads_list.append(F1[4]) 
        for line in hin:
            n +=1
            
            F = line.rstrip('\n').split('\t')
            chr_position = F[0]+':'+F[1]
            
            if chr_position == chr_position_0:
                grouped_reads_list.append(F[4])
                
            else:
                
                consensus, fam_size = consensus_maker(grouped_reads_list, cut_off)
                if read_direction_0 == '-':
                    #reverse comp
                    comp_seq = consensus.translate(trans)
                    revcomp_seq = ''.join(reversed(comp_seq))
                    rec = str(chr_0)+'\t'+ str(position_0)+'\t'+str(sb_length_0)+'\t'+str(read_direction_0)+'\t'+str(fam_size)+'\t'+str(consensus)+'\t'+str(revcomp_seq)+\
                        '\t'+str(sb_direction_0)+'\n'
                    hout.write(rec)
                    
                if read_direction_0 == '+':
                    #
                    rec = str(chr_0)+'\t'+ str(position_0)+'\t'+str(sb_length_0)+'\t'+str(read_direction_0)+'\t'+str(fam_size)+'\t'+str(consensus)+'\t'+str(consensus)+\
                        '\t'+str(sb_direction_0)+'\n'
                    hout.write(rec)
                    
                chr_0 = F[0]
                position_0 = F[1]
                chr_position_0 = F[0]+':'+F[1]
                sb_length_0 = F[2]
                read_direction_0 = F[3]
                sb_direction_0 = F[6]
                grouped_reads_list = [F[4]]
            
            if n == num_lines:
                consensus, fam_size = consensus_maker(grouped_reads_list, cut_off)
                if read_direction_0 == '-':
                    #reverse comp
                    comp_seq = consensus.translate(trans)
                    revcomp_seq = ''.join(reversed(comp_seq))
                    rec = str(chr_0)+'\t'+ str(position_0)+'\t'+str(sb_length_0)+'\t'+str(read_direction_0)+'\t'+str(fam_size)+'\t'+str(consensus)+'\t'+str(revcomp_seq)+\
                        '\t'+str(sb_direction_0)+'\n'
                    hout.write(rec)
                if read_direction_0 == '+':
                    #
                    rec = str(chr_0)+'\t'+ str(position_0)+'\t'+str(sb_length_0)+'\t'+str(read_direction_0)+'\t'+str(fam_size)+'\t'+str(consensus)+'\t'+str(consensus)+\
                        '\t'+str(sb_direction_0)+'\n'
                    hout.write(rec)
    os.remove(tmp_file)
           
                     
if __name__ == "__main__":
	
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output file", required=True)
    parser.add_argument('--cut_off', type=float, default=.8, dest='cut_off', 
						help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
    o = parser.parse_args()
    
    infile = o.infile
    outfile = o.outfile
    cut_off = o.cut_off
    
    run_consensus_maker(infile, outfile, cut_off)


""" 
   
infile='break_BB_341.txt'
outfile='break_cs_BB_341.txt'
cut_off =0.8
"""
