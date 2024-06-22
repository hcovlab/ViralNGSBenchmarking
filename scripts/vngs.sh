#!/bin/bash

Mode="$1"
echo Mode=$Mode
vngs_prefix="$2"
echo vngs_prefix=$vngs_prefix
RawForward="$3"
echo RawForward=$RawForward
RawReverse="$4"
echo RawReverse=$RawReverse
BamReads="$5"
echo BamReads=$BamReads
Adapters="$6"
echo Adapters=$Adapters
AssemblyRef="$7"
echo AssemblyRef=$AssemblyRef

cd /data
/opt/viral-ngs/source/read_utils.py fastq_to_bam --sampleName $vngs_prefix "${RawForward%.*.*}"_r.fastq "${RawReverse%.*.*}"_r.fastq init.bam
/opt/viral-ngs/source/read_utils.py sort_bam init.bam $BamReads coordinate
rm init.bam
/opt/viral-ngs/source/assembly.py assemble_spades $BamReads $Adapters:2:10:7:1:true "${vngs_prefix}".fasta
#/opt/viral-ngs/source/assembly.py assemble_trinity $BamReads $Adapters:2:10:7:1:true "${vngs_prefix}".fasta
if [[ $Mode == "emp_consref" ]]
then
    /opt/viral-ngs/source/assembly.py order_and_orient --min_contig_len=100 --min_pct_contig_aligned=0.3 "${vngs_prefix}".fasta $AssemblyRef "${vngs_prefix}"_ordered.fasta
    /opt/viral-ngs/source/assembly.py impute_from_reference --minLengthFraction=0.3 --minUnambig=0.4 "${vngs_prefix}"_ordered.fasta $AssemblyRef "${vngs_prefix}"_ordered_imputed.fasta --index
    /opt/viral-ngs/source/assembly.py refine_assembly --major_cutoff=0.4 "${vngs_prefix}"_ordered_imputed.fasta $BamReads "${vngs_prefix}"_ordered_imputed_refined.fasta
    if [ ! -s "${vngs_prefix}"_ordered_imputed_refined.fasta ]
    then
        cp "${vngs_prefix}"_ordered_imputed.fasta "${vngs_prefix}"_ordered_imputed_refined.fasta
    fi
else
    /opt/viral-ngs/source/assembly.py order_and_orient "${vngs_prefix}".fasta $AssemblyRef "${vngs_prefix}"_ordered.fasta
    /opt/viral-ngs/source/assembly.py impute_from_reference "${vngs_prefix}"_ordered.fasta $AssemblyRef "${vngs_prefix}"_ordered_imputed.fasta --index
    /opt/viral-ngs/source/assembly.py refine_assembly "${vngs_prefix}"_ordered_imputed.fasta $BamReads "${vngs_prefix}"_ordered_imputed_refined.fasta
fi