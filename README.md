# One-break sleeping beauty
Detection of insertion sites of the Sleeping beauty.
primer-5'SBend- GTGTATGTAAACTTCCGACTTCAACTG+TA

# The place of fastq
../../data/

# Prepare list.txt files.  
1. fastq list   
  fastq_list.txt (commma separate)    
  [0] SeqID in sample_sheet. It is a prefix of fastq.gz file (xxx.fastq.gz, xxx.fastq.gz). Example: 519  
  [1] Tumor number in sample_sheet. It will be a ID for prefix of output fastq name (NX_uniqueID or BB_uniqueID). Example: 59
  [2] NX or BB  
  [3] Group    
Example:  
(Seq ID,	Tumor#,	NXorBB,	Group) the header is not needed.
903,324,NX,35    
904,325,NX,35   
905,326,NX,35   

2. UniqueID and Group list.   
  id_list_all.txt   
  Example  
  110,group  
  111,group
  :

## Mapping
Map the fastq in the fastq folder into the genome.
  
data_dir=../../data  
file=fastq_list.txt  
REFERENCE=../../reference/GRCm38p4_SB.genome.fa
```
run_smap.sh  
```
For single-end, run_smap.sh uses smap.sh.      
   
## STEP2. Breakpoint detection
```
bash run_sb1.sh fastq_list.txt
```
sb1.sh  
(1) breakpoint_detector
```
breakpoint_detectorS.py --infile ${}.bam --outfile ${}_break.txt --re ${re}
```
Use breakpoint_detector.py for paired-end.    
(2) consensus_maker  
Make consensus sequence of reads that have the same insertion position.
```
python consensus_maker.py --infile ./break/${pr}_break.txt --outfile ./break/${pr}_break.cs.txt --cut_off 0.8
```
(3) annotate sb seq and filter the sb inserted positions  
5’- GTGTATGTAAACTTCCGACTTCAACTG ---TA  
1. adjust the insertion position.
2. calculate the sw_match score.
3. calculate the sw_match ratio = sw_match score/soft-clipping size.
4. detect the insertion site bases
5. detect sb motif at break position.
6. extract reads by family size>=3, length of soft-clipping = 22~30 and smith waterman sw_match ratio>0.9.
```
python annot_sbseq.py --input_file ${pr}_break.cs.txt --output_file ${pr}_break.cs.ann.fil.strict.txt
```
## STEP3. Combine files
```
run_sb2.sh
```
Use sb2.sh.
Input file should be "mapping_list.txt".  
The sb2.sh contains 3 steps.  
(1) unify ids   
```
python sb_proc_list.py
```
(2) combine BB and NX 
``` 
python sb_combBBNX_break.py --id {pre} --dir_path ./break
```
Output  
BBNX_002_break.cs.ann.fil.strict.txt   
[0]chr	  
[1]position	  
[2]sb_direction	  
[3]sb_genome_seq  
[4]BBNX_002  

(3) add gene  
``` 
python sb_add_gene.py --input_file ./break/BBNX_${pr}_break.cs.ann.fil.strict.txt --output_file ./break/BBNX_${pr}_break.cs.ann.fil.strict.addgene.txt --gene_file ../reference/gene.bed.gz
```

## STEP3. Make a list for merge.  
```
run_sb3.sh
```
python sb_proc_list2.py  
input_file=./sample/ID_list_all.txt  
ID,Group  
110,1  
111,2  
:  
output_file=ID_list_{Group}.txt

## STEP4. Merge by group, calculate p-value and annotating
Set group list in the run_sb4.sh script.
```
run_sb4.sh 
```
(1)sb_merge_samples.py   
(2)sb_cis.py  
(3)add_mm10repeats.py  
(4)sb_postproc_backTo1.py  
Output:   
..ins1base.humanGene.colrecGene.txt    
new columns  
insertion_count  
position_count  
position	  
human_homolog	 
Freq_in_tcgaColRec  	
Is_cancer_gene  

## make bed file
bash run_sb5.sh  
```
python sb_makebedfile.py --list id_list.txt --outfile output_file 
```

## output file format
1. chr
2. start position (10k window with 5k sliding. The windows are stitched if the p-value is lower than or equal to 0.001.)
3. end position
4. the counts of the insertion sites in the 10k window. The counts are joined by commma in the stitched windows.
5. p-values for each windows.
6. minimum p-value in the stitched windows.
7. Sample IDs which has the insertion in the windows.
8. simpleRepeat score (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/simpleRepeat.txt.gz, https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema)
9. simpleRepeat sequence 
10. Segmental Duplication	(http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=varRep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_doSchema=describe+table+schema)
11. The count of insertions in the stitched window.
12. The total number of the insertion positions.	
13. gene annotations for each insertion position (gene name | direction of the transcript). 
14. The relationship of the direction of Sleeping beauty and the transcript.	
15. human_homolog	
16. Frequency in TCGA Colorectal cancers.	(https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018)
17. Is_cancer_gene

