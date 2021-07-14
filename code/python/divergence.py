## divergence.py
## Ford Fishman
## 6/22/21

import os
from subprocess import check_output, call
from Bio import Seq, SeqIO
from Bio.Align.Applications import MuscleCommandline
import multiprocessing as mp
import numpy as np

ssu_path = "/geode2/home/u030/ffishman/Carbonate/GitHub/MicroSpeciation/data/16S/"
# fasta = ssu_path + 'allseqs_rdp.fasta'
fasta = ssu_path + 'allseqs_silva.fasta'
acc_path = ssu_path + 'blast/allseqs_silva.out'
out_path = ssu_path + 'divergence_silva.csv'
muscle_path = check_output(['which', 'muscle']).decode('utf-8').strip()

seqs =  SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))
# queries = SeqIO.to_dict(SeqIO.parse(fasta, 'fasta'))

info = dict()

def divergence(acc, rdp_id, hits):
    
    seq1 = hits[acc].upper()
    seq2 = hits[rdp_id].upper()

    n = len(seq1)
    assert n == len(seq2), "Unequal sequences: %s, %s" % (acc, rdp_id)

    nucl = set(['A','C','T','G'])

    mismatches = 0

    l = 0 # length where both sequences have nucleotides

    for i in range(n):

        s1i = seq1[i]; s2i = seq2[i]

        if s1i in nucl and s2i in nucl:

            l += 1

            if s1i != s2i: 

                mismatches += 1
    
    div = mismatches/l 

    if l < 400: 
        div = np.nan

    return div

def run_align(seqs:dict, in_path:str, out_path:str, muscle_path=muscle_path):

    with open(in_path, 'w') as f:
        SeqIO.write(seqs, f, 'fasta')

    muscle_cline = MuscleCommandline(muscle_path, input=in_path, out=out_path)
    
    muscle_cline()

    alignment =  SeqIO.to_dict(SeqIO.parse(out_path, 'fasta'))

    return alignment

def map_proc(map_item):
    query = map_item[0]
    hits = map_item[1]
    identities = np.array(map_item[2], dtype=float)
    
    if (identities<97).any():
        return None

    in_path = ssu_path + "temp/" + query + ".fa"
    out_path = ssu_path + "temp/" + query + ".afa"

    alignment = run_align(hits, in_path, out_path)

    divs = list() # divergences

    for ident in alignment.keys():

        if ident != query:

            div = divergence(query, ident, alignment)

            divs.append(div)

    return [query, np.mean(divs)]

query_0 = ''
hits = list()
identities = list()
map_list = list() # dict of lists

with open(acc_path) as f:

    for line in f:

        # parse info
        line = line.strip()
        parts = line.split()
        query = parts[0]
        hit = parts[1]
        ident = parts[3]

        hit_seq = seqs[hit]
        hit_seq.seq = hit_seq.seq.back_transcribe()

        if query == 'conseq':
            continue

        elif query_0 == '':

            query_0 = query
            hits += [seqs[query], hit_seq]
            identities += [ident]

        elif query != query_0:
            map_list.append([query_0, hits, identities])
            # seq_dict[query] = hits
            query_0 = query
        
            hits = [seqs[query], hit_seq]
            identities = [ident]
        
        else:

            hits.append(hit_seq)
            identities.append(ident)

        # if query in info:

        #     info[query] = ['%s\t%s\t%s' % (query, hit, div)]
        
        # else:

        #     info[query].append('%s\t%s\t%s' % (query, hit, div))

pool = mp.Pool(mp.cpu_count())

# print(len(map_list))

print("%s Processors\n" % mp.cpu_count() )

divs = pool.map(map_proc, map_list)

pool.close()

print(divs)

with open(out_path, 'w') as f:

    for entry in divs:

        if not entry is None:

            f.write('%s,%s\n' % (entry[0], entry[1]))
