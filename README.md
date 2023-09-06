# Sleeping beauty Insertion Site Detection
Detection of insertion sites of the Sleeping beauty.  
The sequences around the insertion site should be 
primer-5'SBend- GTGTATGTAAACTTCCGACTTCAACTG+TA - genome

## Input data
Sequencing fastq files are located in the directory ${DATA_DIR}.  
Two librares are prepared using two restriction enzymes, NX or BB for a sample. 

## Preparation of list files.  
### 1. The list of fastq files (fastq_list.txt).   
The file should be a comma-separated file with the following columns:   
- Column[0]: SeqID to identify the fastq file. Use a prefix of the fastq.gz file (e.g., xxx.fastq.gz).
- Column[1]: An unique sample ID to identify the sample. It will be the prefix of the output BAM file (e.g., NX_uniqueID.bam or BB_uniqueID.bam).
- Column[2]: NX or BB
- Column[3]: Group  
Example:   
903,324,NX,35    
904,325,NX,35   
905,326,NX,35  
The header is not needed.   

### 2. The list of Unique IDs and Group names (id_list_all.txt)     
  - Column[0]:  unique id for a sample.
  - Column[1]: group name.  
  Example:   
  110,groupA  
  111,groupA  
  :

## STEP1: Mapping
Map the fastq in ${DATA_DIR} folder to the genome.
```  
data_dir=${DATA_DIR}    
file=${fastq_list.txt}    
REFERENCE=GRCm38p4_SB.genome.fa

run_smap.sh  
```
The GRCm38p4_SB.genome.fa file includes the mouse genome and the Sleeping beauty transposon sequence.  
For single-end sequencing data, the run_smap.sh script utilizes the smap.sh script.      
   
## STEP2: Breakpoint detection
```
bash run_sb1.sh fastq_list.txt
```
This step sb1.sh involves the following three scripts.  
  
(1) Breakpoint detection  
- Detection of breakpoint.  
```
breakpoint_detectorS.py --infile ${pr}.bam --outfile ${pr}_break.txt --re ${re}
```
For paired-end, use breakpoint_detector.py.  

(2) Consensus_maker  
- This step generates consensus sequences from reads with the same insertion position.
```
python consensus_maker.py --infile ./break/${pr}_break.txt --outfile ./break/${pr}_break.cs.txt --cut_off 0.8
```
(3) Annotattion of SB sequences and Filtering the SB inserted positions.  
-  this step performs the following processes:  
3.1. adjust the insertion position.  
3.2. calculate the sw_match score.  
3.3. calculate the sw_match ratio = sw_match score/soft-clipping size.  
3.4. detect the insertion site bases.  
3.5. detect sb motif at break position.  
3.6. extract reads by family size>=3, length of soft-clipping = 22~30 and smith waterman sw_match ratio>0.9.
```
python annot_sbseq.py --input_file ${pr}_break.cs.txt --output_file ${pr}_break.cs.ann.fil.strict.txt
```

## STEP3: Combine files
```
run_sb2.sh
```
In this step, the sb2.sh script is used with the fastq_list.txt.  
This step contains 3 steps.  
(1) Unify IDs   
```
python sb_proc_list.py
```
(2) Combine the result of BB and NX.    
``` 
python sb_combBBNX_break.py --id {pre} --dir_path ./break
```  

(3) Add gene information  
``` 
python sb_add_gene.py --input_file ./break/BBNX_${pr}_break.cs.ann.fil.strict.txt --output_file ./break/BBNX_${pr}_break.cs.ann.fil.strict.addgene.txt --gene_file ../reference/gene.bed.gz
```
Output file format:    
The file prefix starts with "BBNX_" (example, BBNX_002_break.cs.ann.fil.strict.txt ).   
- column[0] chr	  
- column[1] position	  
- column[2] sb_direction	  
- column[3] sb-genome breakpoint sequence  
- column[4] sample name (BBNX_002)  
- column[5] gene  

## STEP4: Create a list used for the merge process.  
```
run_sb3.sh
```
In this step, the following script run with input_file (./sample/ID_list_all.txt).  
The output file will have a name like a id_list_{Group}.txt.
```
python sb_proc_list2.py  
```

## STEP5: Merge files by group, calculate p-value and annotate the gene informations.
Set id_group list in the run_sb4.sh script.
```
run_sb4.sh 
```
- This script includes the following scripts.  
4.1. sb_merge_samples.py   
4.2. sb_cis.py  
4.3. add_mm10repeats.py  
4.4. sb_postproc_backTo1.py  
Output file format:   
The file name will be finish "..ins1base.humanGene.colrecGene.txt"    
 
## STEP6: Make bed file of the insertion postions
```
bash run_sb5.sh  
```
This step involves the following script.  
```
python sb_makebedfile.py --list id_list_{Group}.txt --outfile output_file 
```

## Output file format
1. chr
2. Start position (10k window with 5k sliding. The windows are stitched if the p-value is lower than or equal to 0.001.)
3. End position
4. The counts of the insertion sites in the 10k window. The counts are joined by commma in the stitched windows.
5. p-values for each windows.
6. Minimum p-value in the stitched windows.
7. Sample IDs which has the insertion in the windows.
8. SimpleRepeat score (http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/simpleRepeat.txt.gz, https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=rep&hgta_track=simpleRepeat&hgta_table=simpleRepeat&hgta_doSchema=describe+table+schema)
9. SimpleRepeat sequence 
10. Segmental Duplication	(http://genome.ucsc.edu/cgi-bin/hgTables?db=mm10&hgta_group=varRep&hgta_track=genomicSuperDups&hgta_table=genomicSuperDups&hgta_doSchema=describe+table+schema)
11. The count of insertions in the stitched window.
12. The total number of the insertion positions.	
13. Gene annotations for each insertion position (gene name | direction of the transcript). 
14. The relationship of the direction of Sleeping beauty and the transcript.	
15. Human_homolog	
16. Frequency in TCGA Colorectal cancers.	(https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018)
17. Is_cancer_gene (True when column 16 is present.)

