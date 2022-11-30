#!/bin/bash
#$ -S /bin/bash
#$ -cwd -V
#$ -l mem_req=8G,s_vmem=8G
#$ -pe def_slot 2

REFERENCE=../../reference/GRCm38p4_SB.genome.fa

input=$1
prefix=$2
echo ${input}

gzip -dc ${input} > ${input}.fastq
bwa mem -T 0 ${REFERENCE} ${input}.fastq > ${prefix}.sam
samtools view -Sb ${prefix}.sam > ${prefix}.unsorted.bam
samtools sort ${prefix}.unsorted.bam -o ${prefix}.bam
samtools index ${prefix}.bam
rm ${input}.fastq
rm ${prefix}.sam
rm ${prefix}.unsorted.bam

