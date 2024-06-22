import sys

def fasta_reader(fasta_as_string):
    fasta = []
    for line in fasta_as_string.splitlines():
        if line.startswith('>'):
            header = line.strip()
            fasta.append({header:''})
        elif line.strip():
            fasta[len(fasta)-1][header] += line.strip()
    return fasta

def reverse_complement(sequence: str) -> str:
    complementerBases={'A':'T','C':'G','G':'C','T':'A','-':'-'}
    sequence_rc=''
    for letter in sequence:
        sequence_rc+=(complementerBases[letter])
    return sequence_rc[::-1]

def hamming(text1,text2):
    assert len(text1)==len(text2), 'Different length inputs'
    return sum([1 for idx, elem in enumerate(text1) if elem != text2[idx]])

def approximate_pattern_matching(pattern: str, genome:str,d:int) -> list:
    positions=[]
    idx=0
    while idx+len(pattern)<=len(genome):
        if hamming(genome[idx:(idx+len(pattern))],pattern)<=d :
            positions.append(idx)
        idx+=1
    return positions

def primer_positions(reference:list,primers:list):
    for idx, primdict in enumerate(primers):
        if list(primdict.keys())[0].endswith("F"):
            # Find best position (if multiple select the first)
            d=0
            possible_positions=[]
            while len(possible_positions)==0:
                possible_positions=approximate_pattern_matching(pattern=list(primdict.values())[0],genome=list(reference[0].values())[0],d=d)
                d+=1
            
            # Save positional arguments
            primdict['startpos']=int(possible_positions[0])
            primdict['endpos']=int(possible_positions[0])+len(list(primdict.values())[0])-1
            primdict['strand']="forward"
            if list(primdict.keys())[0].endswith("1F"):
                primdict['nest']="outer"
            else:
                primdict['nest']="inner"
        if list(primdict.keys())[0].endswith("R"):
            complementpattern=reverse_complement(list(primdict.values())[0])
            
            # Find best position (if multiple select last)
            d=0
            possible_positions=[]
            while len(possible_positions)==0:
                possible_positions=approximate_pattern_matching(pattern=complementpattern,genome=list(reference[0].values())[0],d=d)
                d+=1
            
            # Save positional arguments
            primdict['startpos']=int(possible_positions[-1])
            primdict['endpos']=int(possible_positions[-1])+len(list(primdict.values())[0])-1
            primdict['strand']="reverse"
            if list(primdict.keys())[0].endswith("1R"):
                primdict['nest']="outer"
            else:
                primdict['nest']="inner"
    return primers

def bed_from_primers(primers:dict,refname:str):
    with open("Primers.bed", 'w') as f: 
        for idx, primdict in enumerate(primers):
            f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(refname), str(primdict['startpos']), str(primdict['endpos']), str(list(primdict.keys())[0][1:]),str('1' if primdict['nest']=="outer" else '2'),str('+' if primdict['strand']=="forward" else '-')))
            
def tsv_from_primers(primers:dict):
    with open("Primers.tsv", 'w') as f: 
        f.write('%s\t%s\n' % ("name","seq"))
        for idx, primdict in enumerate(primers):
            f.write('%s\t%s\n' % (str(list(primdict.keys())[0][1:]),str(list(primdict.values())[0])))
                 
def insertbed_from_primers(primers:dict,refname:str):
    with open("Inserts.bed", 'w') as f: 
        counter=0
        for idx, primdict in enumerate(primers):
            if idx % 2 == 0:
                insertstart=primdict['endpos']
            if idx % 2 == 1:
                counter+=1
                insertend=primdict['startpos']
                f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(refname), str(insertstart), str(insertend), str(counter),str('1' if primdict['nest']=="outer" else '2'),'+'))        

#### Main
with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    ref_HXB2=fasta_reader(fn.read())
with open(sys.argv[2],'rt',encoding='utf-8-sig') as fn:
    pri_eva=fasta_reader(fn.read())

primers_dict=primer_positions(reference=ref_HXB2,primers=pri_eva)
bed_from_primers(primers=primers_dict,refname=list(ref_HXB2[0].keys())[0].split(" ")[0][1:])
tsv_from_primers(primers=primers_dict)
insertbed_from_primers(primers=primers_dict,refname=list(ref_HXB2[0].keys())[0].split(" ")[0][1:])