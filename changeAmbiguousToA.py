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

def sample_bases(sequence:str):
    changed_genome=""
    for base in list(sequence[0].values())[0]:
        if base in 'RYBDKMHVSWI':
            print("I found a strange character (",base,"), and I will change it to A")
            changed_genome+='A'
        elif base in '()':
            changed_genome+=''
        else:
            changed_genome+=base
    return changed_genome

###########################################x

with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    masterseq=fasta_reader(fn.read())


genome_mod=sample_bases(masterseq)

ofile = open('genome_mod.fasta', "w")  
ofile.write(list(masterseq[0].keys())[0] + "\n" + genome_mod + "\n")
ofile.close()