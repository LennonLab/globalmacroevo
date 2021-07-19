## tree_fasta.py
## Ford Fishman
## 7/15/21

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-wd',  type=str, help="Path to directory with input files")
parser.set_defaults(single=True)
arguments = parser.parse_args()
wd = arguments.wd

acc_file = wd + "/tree_acc.txt"
fasta = wd + "/allseqs.fasta"
output = wd + "/treeseqs.fasta"

seqs =  SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

finalseqs = list()

with open(acc_file,'r') as r:

    for line in r:

        acc = line.strip()
        finalseqs.append(seqs[acc])

with open(output, 'w') as f:
    SeqIO.write(finalseqs, f, 'fasta')