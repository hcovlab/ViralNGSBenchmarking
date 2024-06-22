import sys
import numpy as np
import pandas as pd

def fasta_reader(fasta_as_string):
    fasta = []
    for line in fasta_as_string.splitlines():
        if line.startswith('>'):
            fasta.append('')
        elif line.strip():
            fasta[len(fasta)-1] += line.strip()
    return fasta

                
def compBamLiftover(fastalist,bamframe):
    ma_counter = 0 # multiple alignment counter
    ass_counters = [0]*5 # assembly position counters
    correct, major, ltr, unmapped, extramapped = [0]*5, [0]*5, [0]*5, [0]*5, [0]*5 # assembly position counters
    for read_ind in bamframe.index: # go through every read
        # calculate assembly-specific positions of the mapping
        while ass_counters[0] < bamframe['golden_pos'][read_ind]: # until you reach the position of the golden mapping proceed in sequences
            for ass_ind, assembly in enumerate(fastalist): # go through assemblies
                if assembly[ma_counter] in 'actgnACTGN': # if standard base
                    ass_counters[ass_ind]+=1 # go forward in assembly
                elif assembly[ma_counter] != "-": # if not gap
                    ass_counters[ass_ind]+=1 # go forward in assembly
                    print("Strange character (",assembly[ma_counter],") found in assembly file",assembly[ma_counter])
            ma_counter+=1 # if correct base or gap go forward in aligned sequence
        # compare counters with mapping positions
        for ass_ind, ass_counter in enumerate(ass_counters):
            if ass_counter!=bamframe.iloc[read_ind,ass_ind+1]:
                error_distance=abs(ass_counter - bamframe.iloc[read_ind,ass_ind+1])
                if pd.isna(bamframe.iloc[read_ind,ass_ind+1]):
                    unmapped[ass_ind]+=1
                elif not pd.isna(bamframe.iloc[read_ind,1]) and error_distance<=100:
                    correct[ass_ind]+=1 
                elif not pd.isna(bamframe.iloc[read_ind,1]) and error_distance > 100 and error_distance < 8500:
                    major[ass_ind]+=1
                    print("Found a major (",bamframe.iloc[read_ind,0],") for",bamframe.columns[ass_ind+1],"! The read is mapped to",ass_counter,"instead of",bamframe.iloc[read_ind,ass_ind+1])
                elif not pd.isna(bamframe.iloc[read_ind,1]) and error_distance >= 8500:
                    ltr[ass_ind]+=1
                else:
                    extramapped[ass_ind]+=1
            else:
                correct[ass_ind]+=1
        
    # generate results
    compbam_results = pd.DataFrame(
    {'pos_correct': correct,
     'poserr_major': major,
     'poserr_ltr': ltr,
     'poserr_unmapped': unmapped,
     'poserr_extramapped': extramapped 
    },index=["golden","shiver","smaltalign","viralngs","vpipe"])
    # save results
    with open('MappingAnalysis.txt', 'a') as f:
        compbam_results_string = compbam_results.to_string()
        f.write(compbam_results_string)

########################################################################x

# read input files
with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    fasta_list=fasta_reader(fn.read())

compbams=pd.read_csv(sys.argv[2], sep='\t')
compbams2=compbams.drop(compbams.columns[1], axis=1)
compbams2.columns=["readname","golden","shiver","smaltalign","viralngs","vpipe"]

col_orig=["golden","shiver","smaltalign","viralngs","vpipe"]
colnames=['golden_pos','shiver_pos','smaltalign_pos','viralngs_pos','vpipe_pos']
pattern_start=[':',':',':',':',':']
pattern_end=['=','=','=','=','=']
for colname_ind,colname in enumerate(colnames):
    compbams2[colname]=compbams2[col_orig[colname_ind]].apply(lambda st: st[st.find(pattern_start[colname_ind])+1:st.find(pattern_end[colname_ind])])
    compbams2[colname]=pd.to_numeric(compbams2[colname],errors='coerce')#.astype('Int64')

compbams3=compbams2[["readname"]+colnames].sort_values(by="golden_pos",ignore_index=True)

compBamLiftover(fasta_list,compbams3)

#print(head(compbams3))