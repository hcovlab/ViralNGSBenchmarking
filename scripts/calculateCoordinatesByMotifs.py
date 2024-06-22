import sys
from Bio.Align import PairwiseAligner

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

def calculate_coordinates(motif_pairs:list,sequences:list) -> list:
    aligner = PairwiseAligner(match_score=2, mismatch_score=-2,open_gap_score=-3,extend_gap_score=-1,target_end_gap_score = 0,query_end_gap_score = 0 )
    sequences_cut=[]
    for s_id,sequence in enumerate(sequences):
            # align starter and end motifs
            alignment_start = aligner.align(list(motif_pairs[0].values())[0],list(sequence.values())[0].upper())[0]
            alignment_end = aligner.align(list(motif_pairs[1].values())[0],list(sequence.values())[0].upper())[0]
            # find starting motif first position
            start_id=-1
            start_hit=False
            while not start_hit:
                start_id+=1
                if alignment_start[0][start_id]!='-':
                    start_hit=True
            #find ending motif first position 
            endstart_id=start_id
            endstart_hit=False
            while not endstart_hit:
                endstart_id+=1
                if alignment_end[0][endstart_id]!='-':
                    endstart_hit=True
            #than ending motif last position
            endend_id=endstart_id
            endend_hit=False
            while not endend_hit:
                endend_id+=1
                if alignment_end[0][endend_id:] == len(alignment_end[0][endend_id:]) * alignment_end[0][len(alignment_end[0][endend_id:])-1]:
                    endend_hit=True
            # print sequence length as a proof of the results
            #print("Sequence length:",endend_id-start_id)
            #cut sequence
            if sys.argv[2] == "HXB2_PR_boundaries.fasta":
                coordinates=str(start_id+1)+'-'+str(start_id+296+1)
            else:
                coordinates=str(start_id+1)+'-'+str(endend_id)
    return(coordinates)

###############################################

with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    aligned_sequences=fasta_reader(fn.read())
with open(sys.argv[2],'rt',encoding='utf-8-sig') as fn:
    border_motifs=fasta_reader(fn.read())

genomes_nogaps=delete_gaps(aligned_sequences)
coordinates_range=calculate_coordinates(motif_pairs=border_motifs,sequences=genomes_nogaps)


ofile = open('coordinates.txt', "w")   
ofile.write(coordinates_range)
ofile.close()

