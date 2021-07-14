## taxonomy.py
## Ford Fishman
## 7/13/21

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-wd',  type=str, help="Path to directory with input files")
parser.add_argument('-o', '--out', type=str, help='Output file name')
parser.set_defaults(single=True)
arguments = parser.parse_args()
wd = arguments.wd
out = arguments.out 

def paste_dir(*args):
    return '/'.join(args)

sub = paste_dir(wd, 'subrates.csv') 
blast = paste_dir(wd,'blast/allseqs_silva.out')
fasta = paste_dir(wd,'blast','silva','silva_parsed.fasta')
outfile = paste_dir(wd, out)

seqs =  SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

accs = list()
hits = list()

with open(blast, 'r') as r:

    for line in r:
        line = line.strip()
        parts = line.split('\t')
        query = parts[0]
        hit = parts[1]

        if not query in accs:

            accs.append(query)
            hits.append(hit)

tax_info = list()

for hit in hits:

    record = seqs[hit]
    description = record.description

    for i in range(len(description)):

        char = description[i]

        if char == ' ':
            
            tax_info.append(description[i+1:len(description)])

            break

with open(outfile, 'w') as f:

    for i in range(len(accs)):

        f.write('%s,%s,' % (accs[i], hits[i]) )
        tax = tax_info[i].split(';')
        for ident in tax:

            f.write('%s,' % ident)
        
        f.write('\n')
