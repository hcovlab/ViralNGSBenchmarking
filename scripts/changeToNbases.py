import sys
import random

def fasta_reader(fasta_as_string):
    """
    FASTA reader
    :param: fasta_as_string: fasta formatted string
    :return: FASTA based dictionary with header --> sequence mapping
    """
    fasta = []
    for line in fasta_as_string.splitlines():
        if line.startswith('>'):
            header = line.strip()
            fasta.append({header:''})
        elif line.strip():
            fasta[len(fasta)-1][header] += line.strip()
    return fasta

def ambiguous_bases(sequence:str):
    changed_genome=""
    for base in list(sequence[0].values())[0]:
        if base not in 'actgnACTGN':
            changed_genome+='N'
        else:
            changed_genome+=base
    return changed_genome

###########################################x

with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    masterseq=fasta_reader(fn.read())


genome_mod=ambiguous_bases(masterseq)

ofile = open('genome_mod.fasta', "w")  
ofile.write(list(masterseq[0].keys())[0] + "\n" + genome_mod + "\n")
ofile.close()