## silva_deflines.py
## Ford Fishman
## 7/6/21

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', "--fasta", type=str, help="Path to input fasta")
parser.add_argument('-a', '--acc', type=str, help='Path to accession list')
parser.add_argument('-o', "--out", type=str, help="Path to output fasta")
parser.set_defaults(single=True)
arguments = parser.parse_args()
fasta = arguments.fasta
acc = arguments.acc
output = arguments.out

accs = set()

with open(acc, 'r') as r:

    for line in r:
        accs.add(line.strip())

shouldWrite = True

x = 0

with open(fasta, 'r') as r:

    with open(output, 'w') as f:

        for line in r:

            line = line.strip()

            if line.startswith('>'):
                
                shouldWrite = True
                parts = line.split(' ')
                ident = parts[0].split('.')[0]

                if ident[1:len(ident)] in accs: # remove > character

                    x +=1
               
                    shouldWrite = False
            

            if shouldWrite:
            
                f.write('%s\n' % line)

# print(x)