#!/usr/bin/env python
"""
Breakpoint Detection
By Naoko Iida
2020.03.26

Written for Python 3
Required modules: Pysam, Samtools

Inputs: 
	A position-sorted paired-end BAM file containing reads.

Outputs:
	1. BAM file containing SB sequence
    2. output the chr, break position, soft-clipping and genome(+15 bases) around the break position.

The program starts at the position of the read1  

In the future, this program may be able to autodetect read length.  

usage: breakpoint_detector.py [-h] [--infile INFILE] 
						[--outfile OUTFILE]
						
optional arguments:
	-h, --help            show this help message and exit
	--infile INFILE       input BAM file
	--outfile OUTFILE     output file

"""
import pysam
from argparse import ArgumentParser

def breakpoint_detector(infile, outfile, re):

    #reverse complement
    trans = str.maketrans('ATGCatgc', 'TACGTAGC')
    #infile='BB_1_S1_R1_001.bam'
    #outfile='result.txt'
    # Open the input BAM file
    in_bam_file = pysam.AlignmentFile(f'./bam/{infile}', "rb") 
    # Initialize the iterator
    bam_entry = in_bam_file.fetch(until_eof=True)  
    
    with open(f'./break/{outfile}', 'w') as hout:
    # Start going through the input BAM file, one position at a time.
        for read in bam_entry:
                
            sb_direction = ''
            #print(read)
            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]        
            # skip if not aligned
            if flags[2] == "1": continue
            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue
            
            #get mapq
            mapq = read.mapping_quality
            if mapq <= 10: continue
            if read.is_read1: #<- for paired-end
                #get reference id
                refid = read.reference_name
    
                #get read direction
                if read.is_reverse:
                    direction = '-'
                    tupple = read.cigar[-1]
                else:
                    direction = '+'
                    tupple = read.cigar[0]
                #S 4 soft clipping (clipped sequences present in SEQ)
                if direction == '-' and refid.startswith("chr"): #<----+sb
                    #select soft_clip = True
                    if tupple[0] == 4: 
                        soft_size = tupple[1]
                        aln_length=read.query_alignment_length
                        #reference_end point to one past the last residue. So, break_position is
                        break_posi = read.reference_end-1 
                        a_end = read.query_alignment_end
                        #the size of genome_seq ~15
                        r1 = min([aln_length,15])
                        read.query_sequence[0]
                        #get genome+full of soft-clipping seq. genome_seq+'|'+sb_seq
                        raw_seq = read.query_sequence[a_end-r1:a_end+soft_size]
                        #reverse comp
                        comp_seq = raw_seq.translate(trans)
                        revcomp_seq = ''.join(reversed(comp_seq))
                        #try:
                        #    SA = read.get_tag('SA')
                        #except KeyError:
                        #    continue
                        if re == 'BB':
                            sb_direction = 'sb(->)'
                        elif re == 'NX':
                            sb_direction = 'sb(<-)'
                            
                        rec = str(refid)+'\t'+ str(break_posi)+'\t'+str(soft_size)+'\t'+str(direction)+'\t'+str(raw_seq)+'\t'+str(revcomp_seq)+\
                        '\t'+str(sb_direction)+'\n'
                    else: continue
                
                if direction == '+' and refid.startswith("chr"): #sb+---->
                    if tupple[0] == 4: 
                        soft_size = tupple[1]
                        aln_length=read.query_alignment_length
                        break_posi = read.reference_start
                        a_start = read.query_alignment_start #count 1-start. soft_size
                        #the size of genome_seq ~15
                        r2 = min([aln_length,16])
                        #full of soft-clipping seq + genome. sb_seq+'|'+genome_seq
                        raw_seq = read.query_sequence[0:a_start+r2]
                        #try:
                        #    SA = read.get_tag('SA')
                        #except KeyError:
                        #    continue
                        if re == 'BB':
                            sb_direction = 'sb(<-)'
                        elif re == 'NX':
                            sb_direction = 'sb(->)'
                            
                        rec = str(refid)+'\t'+str(break_posi)+'\t'+str(soft_size)+'\t'+str(direction)+'\t'+str(raw_seq)+'\t'+str(raw_seq)+\
                        '\t'+str(sb_direction)+'\n'
                        
                    else: continue
                hout.write(rec)
    
    	# Close BAM files
        in_bam_file.close()
    	

if __name__ == "__main__":
	
    # Parameters to be input.
    parser = ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile", action="store", dest="outfile", help="output file", required=True)
    parser.add_argument("--re", action="store", dest="re", help="restriction enzyme: BB or NX", required=True)
    o = parser.parse_args()
    
    infile = o.infile
    outfile = o.outfile
    re = o.re
    
    breakpoint_detector(infile, outfile, re)


"""    
infile='BB_341.bam'
main

infile='BB_341.bam'
re = 'BB'
in_bam_file = pysam.AlignmentFile(infile, "rb") 
outfile ='break_BB_341.txt'
"""
