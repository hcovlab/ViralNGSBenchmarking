import sys
from Bio.Align.Applications import MafftCommandline
import tempfile

files = sys.argv[1:]
lines = []

for file in files:
    with open(file, 'r') as ifh:
        lines.append(ifh.read())

with tempfile.NamedTemporaryFile(mode='w') as temp:
    temp.write('\n'.join(lines))
    temp.seek(0)
    mafft_cline = MafftCommandline(input=temp.name)
    stdout, stderr = mafft_cline()

with open('multifasta_aligned.fasta', 'w') as f:
    f.write(stdout)
