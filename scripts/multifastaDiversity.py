import hammingdist
import sys
import numpy as np
import scipy

def fasta_reader_upper(fasta_as_string):
    fasta = []
    for line in fasta_as_string.splitlines():
        if line.startswith('>'):
            if len(fasta)>0:
                fasta[len(fasta)-1][header]=fasta[len(fasta)-1][header].upper()
            header = line.strip()
            fasta.append({header:''})
        elif line.strip():
            fasta[len(fasta)-1][header] += line.strip()
    if len(fasta)>0:
        fasta[len(fasta)-1][header]=fasta[len(fasta)-1][header].upper()
    return fasta

#####################################################

# transform fasta to uppercase
with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    fasta_list=fasta_reader_upper(fn.read())
       
ofile = open("fasta_upper.fasta", "w")
for i in range(len(fasta_list)):
    ofile.write(list(fasta_list[i].keys())[0] + "\n" +list(fasta_list[i].values())[0] + "\n")
ofile.close()

# distance matrix from multifasta
data = hammingdist.from_fasta_large("fasta_upper.fasta")
data.dump("diversity_matrix.csv")

# average distance from distance matrix
distmat=np.loadtxt(open("diversity_matrix.csv", "rb"), delimiter=",", skiprows=0)

distsum=0
for i in range(len(distmat)):
    for j in range(len(distmat[0])):
        if i!=j:
            distsum+=distmat[i,j]
avgdist=distsum/(len(distmat)*(len(distmat)-1))            

print(avgdist)
