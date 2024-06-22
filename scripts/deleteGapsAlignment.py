import sys

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

def hamming(text1,text2):
    assert len(text1)==len(text2), 'Different length inputs'
    return sum([1 for idx, elem in enumerate(text1) if elem != text2[idx]])

def delete_gaps(genomes:str):
    newgenomes=[]
    for g_id,genome in enumerate(genomes):
        idx=0
        newheader=list(genome.keys())[0]
        oldseq=list(genome.values())[0]
        newgenome={newheader:''}
        while idx<len(oldseq):
            if oldseq[idx]!='-':
                newgenome[newheader]+=oldseq[idx]
            idx+=1
        newgenomes.append(newgenome)
    return newgenomes

###############################################

with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    aligned_sequences=fasta_reader(fn.read())

genomes_nogaps=delete_gaps(aligned_sequences)



for idx, nogap_genome in enumerate(genomes_nogaps): # go through every genome  
    ofile = open('refs/'+list(nogap_genome.keys())[0][1:]+'nogap.fasta', "w")  
    ofile.write(list(nogap_genome.keys())[0] + "\n" + list(nogap_genome.values())[0] + "\n")
    ofile.close()

