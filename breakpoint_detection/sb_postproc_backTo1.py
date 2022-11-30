#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 14:04:09 2020

@author: naokoiida

SB_postproc_backTo1.py

1. The begining is the list of candidates.
chr[0]	position_start[1]	position_end[2]	count[3]	p-value[4]	min_p-value[5]	sample[6]	simpleRepeat_score[7]	simpleRepeat_seq	segmentalDups[8]
2. Extract the data of merge_col_BBNX_break.cs.ann.fil.strict.win10k.p.stech.Simprep.Segdups.txt
3. Get a position data from each 
output_file
chr[0] position_start[1] position_end[2] p[4] sample[5] ins_count[6] Pos_count[7] position[8] gene[9]

path_to_data_dir
example: /home/naiida/SleepingBeauty/data/DSS-study
"""

def sb_postproc_backTo1(input_file, path_to_data_dir, group, output_file):
    

    with open(input_file, 'r') as hin, open(output_file, 'w') as hout:

        header = hin.readline().rstrip('\n')
        new_header = header + '\tinsertion_count\tposition_count\tposition\tgene\tsb_gene_direction\n'        
        hout.write(new_header)
        for line in hin:
            
            lie = line.rstrip('\n')
            F = line.rstrip('\n').split('\t')
            chr = F[0]
            start_pos = F[1]
            end_pos = F[2]
            samples = F[6].split(',')
            position_list = []
            gene_list = []
            sb_gene_direction_list = []
            for sample in samples:
                id=sample.split("_")[1].removeprefix('T')
                sample_file = f'{path_to_data_dir}/{group}/{id}/{sample}_break.cs.ann.fil.strict.addgene.txt'
                with open(sample_file) as sin:
                    next(sin)
                    for row in sin:
                        G = row.rstrip('\n').split('\t')
                        if chr == G[0] and int(start_pos) <=int(G[1])<= int(end_pos):
                            position_list.append(G[1])
                            gene_list.append(G[6])
                            strand = G[6].split('|')[1]
                            if strand == '+' and G[2] == '+':
                                sb_gene_direction = 'same'
                            elif strand == '+' and G[2] == '-':
                                sb_gene_direction = 'opposite'                                
                            elif strand == '-' and G[2] == '+':
                                sb_gene_direction = 'opposite'
                            elif strand == '-' and G[2] == '-':
                                sb_gene_direction = 'same'
                            else: sb_gene_direction = '---'
                            
                            sb_gene_direction_list.append(sb_gene_direction)    
            ins_count = len(position_list)
            ins_position_count = len(set(position_list))
            ins_position = ';'.join(set(position_list))
            genes = ';'.join(set(gene_list))
            sb_gene_direction_rec = ';'.join(set(sb_gene_direction_list))
            
            hout.write(lie + '\t' + str(ins_count) +'\t' + str(ins_position_count) +'\t' + str(ins_position) +'\t' + str(genes) +'\t' + str(sb_gene_direction_rec) + '\n')       

def sb_postproc_addHumanhomolog(input_file2, path_to_db):
    
    import pathlib
    
    file_path = pathlib.Path(input_file2)
    name = file_path.stem
    output_file2 = "./merge/"+name+'.humanGene.colrecGene.txt'

    db1 = f'{path_to_db}/Biomart_mm10_Gene_humanHomolog.txt'
    db1_dict = {} 
    with open(db1, 'r') as din1:
        next(din1)
        for line in din1:
            F = line.rstrip('\n').split('\t')
            m_gene = F[0]
            db1_dict[m_gene] = F[1]
    
    db2 = f'{path_to_db}/tcga_pan_can_ColRecMutatedGenes.txt'
    db2_dict = {}
    db3_dict = {}
    with open(db2, 'r') as din2:
        next(din2)
        for line in din2:
            F = line.rstrip('\n').split('\t')
            h_gene = F[0]
            db2_dict[h_gene] = F[5]
            db3_dict[h_gene] = F[6]
            
    with open(input_file2, 'r') as hin, open(output_file2, 'w') as hout:

        header = hin.readline().rstrip('\n')
        new_header = header + '\thuman_homolog\tFreq_in_tcgaColRec\tIs_cancer_gene\n'
        hout.write(new_header)
        
        for line in hin:
            lie = line.rstrip('\n')
            R = line.rstrip('\n').split('\t')
            genes = R[13].split(';')
            
            homolog_genes_list = []
            for gene in genes:
                if gene == '---|---': continue
                else:
                    q = gene.split('|')[0]
                    homolog_gene = db1_dict.get(q, '-')
                    homolog_genes_list.append(homolog_gene)
            homolog_genes = ','.join(homolog_genes_list)
            freq_list = []
            cancer_list = []
            for h in homolog_genes.split(','):
                freq_list.append(db2_dict.get(h, '-'))
                cancer_list.append(db3_dict.get(h,'-'))
            freq = ','.join(freq_list)
            cancer = ','.join(cancer_list)
 
            rec = lie +'\t'+ homolog_genes +'\t'+freq+'\t'+cancer + '\n'
            hout.write(rec)   

if __name__== "__main__":
    import argparse
    import pathlib
    
    parser = argparse.ArgumentParser() #make a parser
    
    parser.add_argument("--input_file", action="store", metavar = "input_file", default = None, type = str,
                            help = "a query file.") 
    parser.add_argument("--path_to_data_dir", action="store", metavar = "path_to_data_dir", default = None, type = str,
                            help = '/path/to/folder/ for /group/id/BBNX_id_break.cs.ann.fil.strict.txt') 
    parser.add_argument("--path_to_db",action="store", metavar = "path_to_db", default = '.', type = str,
                            help = "/path/to/database") 
    parser.add_argument("--group",action="store", metavar = "group", default = '.', type = str,
                            help = "group name")             
    args = parser.parse_args()    
    
    input_file = args.input_file
    path_to_data_dir = args.path_to_data_dir
    path_to_db = args.path_to_db
    group = args.group
        
    file_path = pathlib.Path(input_file)
    name = file_path.stem
    output_file = "./merge/" + name+'.ins1base.txt'
    
    sb_postproc_backTo1(input_file, path_to_data_dir, group, output_file)
    sb_postproc_addHumanhomolog(output_file, path_to_db)





