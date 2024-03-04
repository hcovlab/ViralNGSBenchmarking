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

def matchCoordinates(fastalist,vcflist):
    for vcf_index,vcf in enumerate(vcf_list): # go through every vcf file
        ma_counter=0 # multiple alignment counter
        fasta_counter=0 # original fasta sequence counter
        fasta=list(fasta_list[vcf_index].values())[0] # extract fasta sequence
        for index,vcf_row in vcf.iterrows(): # go through every line in vcf
            while fasta_counter<int(vcf_row['POS']): # while there are variants left in the vcf file
                if fasta[ma_counter] in 'actgnACTGN': # if not standard base
                    fasta_counter+=1 # go forward in original fasta sequence
                elif fasta[ma_counter] != "-": # if not gap
                    fasta_counter+=1 # go forward in original fasta sequence
                    print("Strange character (",fasta[ma_counter],") found in VCF file",vcf_index,"in sequence position",vcf_row['POS'])
                ma_counter+=1 # if correct base or gap go forward in aligned sequence
            if fasta[ma_counter-1].upper()==vcf_row['REF'].upper(): # if vcf file position and sequence in consensus match
                vcf_row['POS']=ma_counter # adjust coordinates according to the multiple alignment
            else:
                print("Bases in VCF and consensus do not match for sequence",vcf_index,"in position",vcf_row['POS'])
    return(vcf_list)

def prediction_arrays(fasta_list,vcf_list):
    positions=list(range(len(list(fasta_list[0].values())[0])))
    predlist=[]
    for vcf in vcf_list:
        predlist_vcf=[]
        for position in positions:
            if position+1 in list(vcf['POS']):
                print("possible variant position")
                is_minvar=False
                vcf_positions=np.where( vcf['POS'] == (position+1))[0]
                for vcf_position in vcf_positions:
                    print("Allele frequency:",float(vcf.iloc[vcf_position]['AF']))
                    if float(vcf.iloc[vcf_position]['AF'])!=1.0:
                        is_minvar=True
                if is_minvar:
                    predlist_vcf.append(True)
                    print("variant found")
                else:
                    predlist_vcf.append(False)    
            else:
                predlist_vcf.append(False)
        predlist.append(predlist_vcf)
    return predlist

########################################################################x

# read input files
with open(sys.argv[1],'rt',encoding='utf-8-sig') as fn:
    fasta_list=fasta_reader(fn.read())

vcf_files=sys.argv[2:]
vcf_list=[]
for vcf_file in vcf_files:
    vcf_list.append(VCF.dataframe(vcf_file, large=False))
print(vcf_list)

# match coordinates
vcf_list=matchCoordinates(fasta_list,vcf_list)

# create prediction lists
predlist=prediction_arrays(fasta_list,vcf_list)

# extract sequence names from fasta list
seqs=[]
for fasta in fasta_list[1:]:
    seqs.append(list(fasta.keys())[0])

# calculate metrics
df = pd.DataFrame(columns=['Accuracy','Precision','Sensitivity','Specificity','F1-score'], index=seqs)
for pred_id,pred in enumerate(predlist[1:]):
    df.loc[seqs[pred_id]] = pd.Series({
        'Accuracy':metrics.accuracy_score(y_true=predlist[0],y_pred=pred),
        'Precision':metrics.precision_score(y_true=predlist[0],y_pred=pred,zero_division=0.0),
        'Sensitivity':metrics.recall_score(y_true=predlist[0],y_pred=pred,zero_division=0.0),
        'Specificity':metrics.recall_score(y_true=predlist[0],y_pred=pred,pos_label=0,zero_division=0.0),
        'F1-score':metrics.f1_score(y_true=predlist[0],y_pred=pred,zero_division=0.0)
        })

# save results
with open('VariantAnalysis.txt', 'a') as f:
    df_string = df.to_string()
    f.write(df_string)