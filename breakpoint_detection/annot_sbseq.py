#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 15:09:45 2020

@author: naokoiida

5’- GTGTATGTAAACTTCCGACTTCAACTG ---TA

1. adjust the insertion position.
2. calculate the sw_match score.
3. calculate the sw_match ratio = sw_match score/soft-clipping size.
4. detect the insertion site bases
5. detect sb motif at break position.
6. extract reads by family size>=3, length of soft-clipping = 12~30 and smith waterman sw_match ratio>0.9.

column of input file
            1. chr
            2. position
            3. soft-clipping length
            4. read_direction
            5. reads
            6. seq
            7. sb_direction

smith waterman 
ref ='GTGTATGTAAACTTCCGACTTCAACTGTA'
query = 'TGTATGTAAACTTCCGACTTCAACTCTGTATCAGTGAGTACAT'
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
alignment = sw.align(query, ref)
alignment.dump()

column of output file
            1. chr
            2. position (after adjustment)
            3. original position
            4. family size
            5. sb direction
            6. seq (after adjustment)
            7. soft-clipping length
            8. sw_match_score
            9. sw_match_ratio
            10. the sb seq at the insertion site.
            11. the genome seq at insertion site.
            12. ori_sb|genome_seq

"""


from argparse import ArgumentParser
import swalign

def annot_sbinsert(input_file, output_file):
    #choose your own values here… 2 and -1 are common.
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    
    #5’- GTGTATGTAAACTTCCGACTTCAACTG ---TA
    seq_dict = {'CGACTTCA': -4,'GACTTCAA': -3,'ACTTCAAC': -2,'CTTCAACT': -1,'TTCAACTG': 0,
                'CCGACTTC':-5,'TCCGACTT':-6,'TTCCGACT':-7,'CTTCCGAC':-8,'ACTTCCGA':-9,
                'AACTTCCG':-10,'AAACTTCC':-11}
    
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:
        next(hin)
        header = '\t'.join(['chr','position(0-start)','ori_position','family_size','sb_direction','adj_seq','soft_clip_len','sw_match','sw_match_ratio','sb-seq|','|genome2bases'])
        hout.write(header+'\n')
        
        for line in hin:

            F = line.rstrip('\n').split('\t')
            #seq ='GTGTATGTAAACTTCCGACTTCAACTGTAATTCTCTGAATGG'
            chr = F[0]
            position = F[1]
            read_direction = F[3]
            reads = F[4]
            seq = F[6]
            sb_length = int(F[2])
            sb_direction = F[7]
            break_motif = seq[sb_length-8:sb_length]
            j = 0
            #check sb motif in nearby break position.
            if break_motif in seq_dict.keys():
                sb_motif = '+'
                adj = seq_dict.get(break_motif)
                if adj == -4 and seq[sb_length-8:sb_length+4] == 'CGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+4:sb_length+6]
                elif adj == -3 and seq[sb_length-8:sb_length+3] == 'GACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+3:sb_length+5]
                elif adj == -2 and seq[sb_length-8:sb_length+2] == 'ACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+2:sb_length+4]
                elif adj == -1 and seq[sb_length-8:sb_length+1] == 'CTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+1:sb_length+3]
                elif adj == 0 and seq[sb_length-8:sb_length] == 'TTCAACTG':
                    j = 1
                    genome2base = seq[sb_length:sb_length+2]
                
                #adj from -5 to -11  #seq ='GTGTATGTAAACTTCCGACTTCAACTGTAATTCTCTGAATGG'
                elif adj == -11 and seq[sb_length-8:sb_length+11] == 'AAACTTCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+11:sb_length+13]
                elif adj == -10 and seq[sb_length-8:sb_length+10] == 'AACTTCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+10:sb_length+12]
                elif adj == -9 and seq[sb_length-8:sb_length+9] == 'ACTTCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+9:sb_length+11]
                elif adj == -8 and seq[sb_length-8:sb_length+8] == 'CTTCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+8:sb_length+10]
                elif adj == -7 and seq[sb_length-8:sb_length+7] == 'TTCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+7:sb_length+9]
                elif adj == -6 and seq[sb_length-8:sb_length+6] == 'TCCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+6:sb_length+8]
                elif adj == -5 and seq[sb_length-8:sb_length+5] == 'CCGACTTCAACTG':
                    j = 1
                    genome2base = seq[sb_length+5:sb_length+7]
                    
                if j == 1:
                    #adj_seq ori_sb_genome_seq = seq[0:sb_length] + '|' + seq[sb_length:]
                    adj_seq = seq[0:sb_length+int(adj)*-1] + '|' + seq[sb_length+int(adj)*-1:]
                    if read_direction == '-':
                        adj_position = int(position) + int(adj)
                    if read_direction == '+':
                        adj_position = int(position) + int(adj)*-1
                else:
                    adj_position = int(position)
                    sb_motif = '-'
                    adj_seq = seq[0:sb_length] + '|' + seq[sb_length:]
                    genome2base = seq[sb_length:sb_length+2]
            else:
                #adj_position = original position
                adj_position = int(position)
                sb_motif = '-'
                genome2base = seq[sb_length:sb_length+2]
                adj_seq = seq[0:sb_length] + '|' + seq[sb_length:]
                    
            #swalign
            alignment = sw.align(seq, 'GTGTATGTAAACTTCCGACTTCAACTG')
            sw_match = alignment.matches
            #sw_cigar = alignment.cigar 
            #ratio sw_match score/soft-clipping length
            sw_ratio = round(float(sw_match/sb_length),2)
            
            #'chr','position(0-start)','ori_position','family_size','sb_direction','adj_seq','soft_clip_len','sw_match','sw_match_ratio','sb-seq|','|genome2bases','ori_sb|genome_seq'
            rec = str(chr)+'\t'+str(adj_position)+'\t'+str(position)+'\t'+str(reads)+'\t'+str(sb_direction)+'\t'+str(adj_seq)+'\t'+ \
                    str(sb_length)+'\t'+str(sw_match)+'\t'+str(sw_ratio)+'\t'+ \
                    str(sb_motif)+'\t'+str(genome2base)+'\n'         
                    # filtering: reads>=3, sb_length =22~30 and sw_match_ratio>0.9 PASS
            if int(reads)>=3 and float(sw_ratio)>0.9 and int(sb_length)>=22 and int(sb_length)<=30 and str(sb_motif)=='+' and str(genome2base)=='TA':
                hout.write(rec)
            else: continue        
            
if __name__ == "__main__":
	
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--input_file", action="store", dest="input_file", help="input file", required=True)
    parser.add_argument("--output_file", action="store", dest="output_file", help="output file strict filtering", required=True)
    
    o = parser.parse_args()
    input_file = o.input_file
    output_file = o.output_file
    
    annot_sbinsert(input_file, output_file)

