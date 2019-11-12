#!/bin/bash
#11/11/19
#script to immitate Macio :)
#converts fastq file to fasta

fastq=/home/jason/practice/fastq2fasta/test.fastq

(head -1 $fastq ; grep -A 1 "@" $fastq | sed '/--/d;/@/d') | sed 's/@/>/' > out.fasta
cat out.fasta
