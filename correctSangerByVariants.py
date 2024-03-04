import sys
import VCF
from pysam import VariantFile
import numpy as np
from sklearn import metrics
import pandas as pd

def fasta_reader(fasta_as_string):
    fasta = []
    for line in fasta_as_string.splitlines():
        if line.startswith('>'):
            header = line.strip()
            fasta.append({header:''})
        elif line.strip():
            fasta[len(fasta)-1][header] += line.strip()
    return fasta

def correct_fasta(fastaseq,vcf):
    fastaseq=list(fastaseq[0].values())[0]
    changed_genome=""
    for position, base in enumerate(fastaseq):
        if (position+1 in list(map(int, vcf['POS']))):
            vcf_positions=np.where( vcf['POS'] == str(position+1))
            newbase=''
            refbase=''
            for vcf_position in vcf_positions[0]:
                vcf_position=int(vcf_position)
                refbase=vcf.iloc[vcf_position]['REF']
                if (float(vcf.iloc[vcf_position]['AF'])>0.7 and refbase=="A"):
                    newbase=vcf.iloc[vcf_position]['ALT']
                    print("I changed position",position+1,"in the Sanger sequence from",refbase,"to",newbase)
            if newbase!='':
                changed_genome+=newbase
            else:
                changed_genome+=refbase
        elif base in '()':
            changed_genome+=''
        else:
            changed_genome+=base
    return changed_genome

########################################################################x

# read input files
with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    fastaseq=fasta_reader(fn.read())

vcf_file=sys.argv[2]
vcf_dataframe=VCF.dataframe(vcf_file, large=False)
print(vcf_dataframe)

# correct fasta
fasta_corrected=correct_fasta(fastaseq,vcf_dataframe)

ofile = open('Sanger_consref_corrected.fasta', "w")  
ofile.write(">corrected_Sanger" + "\n" + fasta_corrected + "\n")
ofile.close()