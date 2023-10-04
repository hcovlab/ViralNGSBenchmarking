import textwrap
import sys

def readHMM(filename:str):
    hmmLines=[]
    with open(filename) as f:
        Lines = f.readlines()
        for line in Lines:
            hmmLines.append(line.strip().split())
    return hmmLines
          
def  hmmToFasta(hmmLines):
    fastaseq=""
    linecounter=1
    for hmmLine in hmmLines:
        if hmmLine[0]=="M":
            break
        if linecounter==3:
            basekey=hmmLine
        elif linecounter>3:
            maxInd=-1
            maxVal=-1
            for ind,prop in enumerate(hmmLine):
                if float(prop)>maxVal:
                    maxInd=ind
                    maxVal=float(prop)
            fastaseq=fastaseq+basekey[maxInd]
        linecounter+=1    
    return fastaseq

def fastaWriter(filename,header,seq):
    with open(filename, "w") as outfile:
        outfile.write(header + "\n")
        outfile.write("\n".join(textwrap.wrap(seq, 60)))
        outfile.write("\n")
              
### Main 
print("Writing hmm profile consensus to fasta sequence...")
hmmLines=readHMM(sys.argv[1])
fastaseq=hmmToFasta(hmmLines)
fastaWriter("vpipe_consensus_prefixed.fasta",">currsample-currdate",fastaseq)
print("...and the writing was successful! Have a nice day!")
       