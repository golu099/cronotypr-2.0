# test1 split up fasta file using coordinates
import sys
import Bio
from Bio import SeqIO
#input_stuff
with open(sys.argv[1], 'r') as f:
        coords = {}
        for line in f:
                coords[line.strip().split('\t')[0]] = (int(line.strip().split('\t')[1]), int(line.strip().split('\t')[2]))

#Matrix formed
with open(sys.argv[2], 'r') as fasta:
        with open(sys.argv[3], 'w') as out:
                for line in fasta:
                        if line.startswith('>'):
                            header = line.strip().split('>')[1]
                            out.write('>' + header + '\n' + next(fasta).strip()[coords[header][0]:coords[header[1]]])
