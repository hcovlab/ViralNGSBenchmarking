#!/bin/bash

while read folder read seq
do
    mkdir -p $folder
    cd $folder
    SRR=$(pysradb srx-to-srr $read | grep -o 'SRR[0-9]*')
    fasterq-dump --split-files $SRR -O .
     rename 's/_1/_R1/g' *.fastq; rename 's/_2/_R2/g' *.fastq
    efetch -db nuccore -id $seq -format fasta > ${seq}_consref.fa
    cd ..
done < $1