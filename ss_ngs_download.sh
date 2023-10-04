#!/bin/bash

while read folder read seq
do
    mkdir -p $folder
    cd $folder
    SRR=$(pysradb srx-to-srr $read | grep -o 'SRR[0-9]*')
    fasterq-dump --split-files $SRR -O .
    efetch -db nuccore -id $seq -format fasta > ${seq}_consref.fa
    cd ..
done < $1