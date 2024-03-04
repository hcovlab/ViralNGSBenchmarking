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
    dist=0
    for idx in range(len(text1)):
        if text2[idx]=='N' or text1[idx]=='N':
            dist+=0.3
        elif text1[idx] != text2[idx]:
            dist+=1
    return dist

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

def cut_edges(refseq:list,sequences:list,HXB2_start:int, HXB2_end:int) -> list:
    #convert coordinates to 0-based
    HXB2_start-=1
    HXB2_end-=1
    sequences_cut=[]
    for s_id,sequence in enumerate(sequences):
        # try:
            #print("New seq!")
            # set aligner
            if str(list(sequence.values())[0]).count('N')>50:
                aligner = PairwiseAligner(match_score=2, mismatch_score=-2,open_gap_score=-3,extend_gap_score=-2,target_end_gap_score = 0,query_end_gap_score = 0 )
            else:
                aligner = PairwiseAligner(match_score=2, mismatch_score=-2,open_gap_score=-3,extend_gap_score=-1,target_end_gap_score = 0,query_end_gap_score = 0 )
            # align HXB2 and sequence
            alignment = aligner.align(list(refseq[0].values())[0],list(sequence.values())[0].upper())[0]
            #print(alignment)
            # find start position
            alignment_counter=-1
            ref_id=-1
            start_hit=False
            while not start_hit:
                alignment_counter+=1
                if alignment[0][alignment_counter]!='-':
                    ref_id+=1
                if ref_id==HXB2_start:
                    start_hit=True
            start_id=alignment_counter
            #print("START: HXB2 next five:", list(refseq[0].values())[0][HXB2_start:(HXB2_start+5)], "alignment next five",alignment[0][start_id:(start_id+5)])
            # find end position     
            end_id=start_id
            end_hit=False
            while not end_hit:
                alignment_counter+=1
                if alignment[0][alignment_counter]!='-':
                    ref_id+=1
                if ref_id==HXB2_end:
                    end_hit=True
            end_id=alignment_counter
            #print("END: HXB2 previous five:", list(refseq[0].values())[0][(HXB2_end-5):HXB2_end], "alignment previous five",alignment[0][(end_id-5):end_id])   
            # cut and save
            newseq={}
            newseq[list(sequence.keys())[0]]=alignment[1][start_id:end_id]
            #print("Reference sequence:",list(refseq[0].values())[0][HXB2_start:HXB2_end+1])
            #print("New sequence",alignment[1][start_id:end_id+1])
            sequences_cut.append(newseq)
            #print('Good!')
        # except:
        #     newseq={}
        #     newseq['errorseq']='N'*((HXB2_end+1)-HXB2_start)
        #     sequences_cut.append(newseq)
    return(sequences_cut)

###############################################

with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    aligned_sequences=fasta_reader(fn.read())
with open(sys.argv[2],'rt',encoding='utf-8-sig') as fn:
    HXB2=fasta_reader(fn.read())

genomes_nogaps=delete_gaps(aligned_sequences)
cropped_genomes=cut_edges(refseq=HXB2,sequences=genomes_nogaps,HXB2_start=int(sys.argv[3]),HXB2_end=int(sys.argv[4]))
cropped_genomes=delete_gaps(cropped_genomes)

ofile = open(sys.argv[5], "w")
for idx, cropped_genome in enumerate(cropped_genomes): # go through every genome    
    ofile.write(list(cropped_genome.keys())[0] + "\n" + list(cropped_genome.values())[0] + "\n")
ofile.close()

