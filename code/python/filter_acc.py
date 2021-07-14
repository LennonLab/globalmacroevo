## filter_acc.py
## Ford Fishman
## 6/16/21

ssu_path = "/geode2/home/u030/ffishman/Carbonate/GitHub/MicroSpeciation/data/16S/"
main_path = ssu_path + "allseqs.out"
bac_path = ssu_path + "undersampled_bacteria_acc.txt" # uncultured bacteria
fun_path = ssu_path + "fungi_acc.txt" # fungi
arc_path = ssu_path + "archaea_acc.txt" # archaea

all_accs = list() # accessions

with open(main_path) as r:
    
    for line in r:

        line = line.strip()
        parts = line.split('\t')
        acc = parts[0] # grab accession
        all_accs.append(acc)

def read_acc_file(path):

    with open(path) as r:
        accs = set()
    
        for line in r:

            acc = line.strip()
            accs.add(acc)
    
    return accs

# read in exceptions
bac_accs = read_acc_file(bac_path)
fun_accs = read_acc_file(fun_path)
arc_accs = read_acc_file(arc_path)

final_acc = list()

for acc in all_accs:

    if not acc in bac_accs and not acc in fun_accs and not acc in arc_accs:

        final_acc.append(acc)

# output
stseqs_path = ssu_path + "stseqs_acc.txt" # standard sequences

with open(stseqs_path, "w") as f:

    for acc in final_acc:

        f.write("%s\n" % acc)
