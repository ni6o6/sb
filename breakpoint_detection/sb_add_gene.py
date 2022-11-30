#! /usr/bin/env python

#python add_gene.py inputfile outputfile gene.bed.gz

import pysam


def sb_addgene(input_file, output_file, gene_file):
    
    gene_tb = pysam.TabixFile(gene_file)
    
    with open(input_file, 'r') as hin, open(output_file, 'w') as hout :
    
        header = hin.readline().rstrip('\n')
    
    
        new_header = header + '\tGene\tStrand\n'
        hout.write(new_header)
        
        for line in hin:
            lie = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
    
            tchr = F[0]
            tpos = F[1]
    
            tabixErrorFlag = 0
            try:
                records = gene_tb.fetch(tchr, int(tpos) - 1, int(tpos) + 1)
            except Exception as inst:
                print("%s: %s" % (type(inst), inst.args))
                tabixErrorFlag = 1
    
            gene_list = []
            strand_list = []
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')
                    gene_list.append(record[3])
                    strand = record[3]+'|'+record[5]
                    strand_list.append(strand)
    
            if len(gene_list) == 0: gene_list.append("---")
            if len(strand_list) == 0: strand_list.append("---|---")
            gene = ';'.join(list(set(gene_list)))
            strand = ';'.join(list(set(strand_list)))
    
            rec = lie +'\t'+ gene +'\t'+ strand + '\n'
            hout.write(rec)
        
if __name__== "__main__":
    import argparse
    
    parser = argparse.ArgumentParser() #make a parser
    parser.add_argument("--input_file", action="store", dest="input_file",
                            help = "input file", required=True)
    parser.add_argument("--output_file", action="store", dest="output_file",
                            help = "output file", required=True)
    parser.add_argument("--gene_file", action="store", dest="gene_file",
                            help = "bed file for tabix", required=True)

    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    gene_file = args.gene_file
 
    sb_addgene(input_file, output_file, gene_file)

