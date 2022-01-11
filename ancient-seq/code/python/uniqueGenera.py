import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--file', default='allseqs.out', type=str, help="Desired name of input path")
arguments = parser.parse_args()
filename = arguments.file

species = dict()

with open(filename) as r:
    
    for line in r:

        line = line.strip()
        parts = line.split('\t') 

        query = parts[0]


        if not query in species:

            species[query] = [[parts[6]],[parts[3]]] #[[species],[perc. ident]]

        else:
            species[query][0].append(parts[6]) # species
            species[query][1].append(parts[3]) # perc. ident

final_genera = dict()

for query, parts in species.items():
    spp = parts[0]; ident = parts[1]
    genera = [sp.split(" ")[0] for sp in spp] # take the genus

    counts = dict()

    for i in range(len(genera)):
        genus = genera[i]
        pident = ident[i]

        if genus in counts: # scale results by percent identity
            counts[genus] += pident

        else:
            counts[genus] = pident
    
    # pick the genus with the highest weighted count
    genus = max(counts, key=counts.get)
    final_genera[query] = genus

with open("query_genera.txt", "w") as f:

    for query, genus in final_genera.items():

        f.write("%s\t%s\n" % (query, genus))

unique_genera = set(final_genera.values())

with open("unique_genera.txt", "w") as f:

    for genus in unique_genera:
        f.write("%s\n" % genus)
    