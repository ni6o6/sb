# sleeping beauty insertion
screening of insertion sites of the Sleeping beauty.
primer-5'SBend- GTGTATGTAAACTTCCGACTTCAACTG+TA


## Mapping
Please save your fastq files in the fastq folder.
.sh
```
FOLDER=$1
REFERENCE=./reference/GRCm38p4_SB.genome.fa

mkdir ./${FOLDER}/bam
mkdir ./${FOLDER}/break

for INPUT in `\find ./${FOLDER}/fastq/ -name '*.fastq'`; do

echo ${INPUT}
IFS='/' read -r -a file <<< "${INPUT}"
IFS='_' read -r -a array <<< "${file[3]}"

file_name=${file[3]}
pr=${array[0]}_${array[1]}
echo ${file_name}
echo ${pr}

cd ./${FOLDER}/fastq/
gzip -dc ./${FOLDER}/fastq/${file_name} > ./${FOLDER}/fastq/${pr}.fastq
bwa mem -T 0 ${REFERENCE} ./${FOLDER}/fastq/${pr}.fastq > ./${FOLDER}/bam/${pr}.sam
samtools view -Sb ./${FOLDER}/bam/${pr}.sam > ./${FOLDER}/bam/${pr}.unsorted.bam
samtools sort ./${FOLDER}/bam/${pr}.unsorted.bam -o ./${FOLDER}/bam/${pr}.bam
samtools index ./${FOLDER}/bam/${pr}.bam
rm ./${FOLDER}/fastq/${pr}.fastq
```
paired-end
```
SB_2020n1
FOLDER=$1
REFERENCE=./reference/GRCm38p4_SB.genome.fa  

mkdir ./${FOLDER}/bam
mkdir ./${FOLDER}/break

for INPUT in `\find ./${FOLDER}/fastq/*/ -name '*_1.fq.gz'`; do

echo ${INPUT}  
IFS='/' read -r -a file <<< "${INPUT}"  
IFS='_' read -r -a array <<< "${file[4]}"  

folder_s=${array[0]}_${array[1]}_${array[2]}
fq1=${array[0]}_${array[1]}_${array[2]}_${array[3]}_${array[4]}_1.fq.gz fq2=${array[0]}_${array[1]}_${array[2]}_${array[3]}_${array[4]}_2.fq.gz
pr=${array[0]}_${array[1]}  
echo ${fq1}
echo ${fq2}
echo ${pr}

#cd ./${FOLDER}/fastq/${folder_s}  
gzip -dc ./${FOLDER}/fastq/${folder_s}/${fq1} > ./${FOLDER}/fastq/${folder_s}/${pr}_1.fastq  
gzip -dc ./${FOLDER}/fastq/${folder_s}/${fq2} > ./${FOLDER}/fastq/${folder_s}/${pr}_2.fastq  
bwa mem -T 0 ${REFERENCE} ./${FOLDER}/fastq/${folder_s}/${pr}_1.fastq ./${FOLDER}/fastq/${folder_s}/${pr}_2.fastq > ./${FOLDER}/bam/${pr}.sam  
samtools view -Sbh ./${FOLDER}/bam/${pr}.sam > ./${FOLDER}/bam/${pr}.unsorted.bam  
samtools sort ./${FOLDER}/bam/${pr}.unsorted.bam -o ./${FOLDER}/bam/${pr}.bam  
samtools index ./${FOLDER}/bam/${pr}.bam  
rm ./${FOLDER}/fastq/${folder_s}/${pr}_1.fastq
rm ./${FOLDER}/fastq/${folder_s}/${pr}_2.fastq
```
## breakpoint_detector
```
breakpoint_detector.py --infile ${pr}.bam --outfile ${pr}_break.txt --re BB
```
* breakpoint_detectorS.py for single-end.

input bam is in bam folder.<br>
output resultant file in break folder.

## consensus_maker
make consensus sequence of reads that have the same insertion position.
```
python consensus_maker.py --infile ${pr}_break.txt --outfile ${pr}_break.cs.txt --cut_off 0.8
```

## annotate sb seq and filter the sb inserted positions
5’- GTGTATGTAAACTTCCGACTTCAACTG ---TA

1. adjust the insertion position.
2. calculate the sw_match score.
3. calculate the sw_match ratio = sw_match score/soft-clipping size.
4. detect the insertion site bases
5. detect sb motif at break position.
6. extract reads by family size>=3, length of soft-clipping = 24~30 and smith waterman sw_match ratio>0.9.

```
python annot_sbinsert.py --infile ${pr}_break.cs.txt --out_pre ${pr}_break.cs
```

## combine and merge samples
prepare a list file with sample number.
   example
   110
   111

## combBBNX_break

folder name = ./break/
input file
    BB_{sample}
    NX_{sample}

sh SB_comb_merge.sh list.txt
    python SB_combBBNX_break.py --sample {sample}
    python SB_merge_sample.py --list {list.txt} --merge_name

## Post processing
merge samples

'''
python 
FOLDER=$1 #folder name<br>
SIZE=10000 #10k,100k<br>
WIN=win10bin5k #use for file name<br>
SB_fill.py<br>
SB_cis_slidingwindow.py<br>
SB_calpvalue_win.py<br>
