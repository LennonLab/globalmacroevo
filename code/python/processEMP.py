## Ford Fishman
## processEMP.py
## 5/21/21

import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial

def map_proc(line, site_labels)->bool:
    
    line = line.strip()
    parts = line.split("\t")
    otu = parts[1:-1]

    # vals = np.array([ -1 if env_dict[sample]=="Free-living" else 1 for sample in samples ]) 

    abundance = np.array(otu, dtype='f')
    
    pa = np.heaviside(abundance, 0) # presence absence

    env_type = site_labels * pa # -1 for free-living, 0 for absent, 1 for host-associated

    # booleans
    obligate_ha = (env_type >= 0).all() #obligate host associated
    obligate_fa = (env_type <= 0).all() #obligate free-living
    ha = (env_type > 0).any() # host associated
    fa = (env_type < 0).any() # free-living
    majority_ha = np.mean(env_type) > 0 # majority host-associated
    majority_fa = np.mean(env_type) < 0 # majority free-living
       
    return [ha, majority_ha, obligate_ha, fa, majority_fa, obligate_fa]

emp_path = "/geode2/home/u030/ffishman/Carbonate/GitHub/MicroSpeciation/data/EMP/"
sbs_path = emp_path + "emp_sbs.txt"
meta_path = emp_path + "emp_qiime_mapping_release1_20170912.tsv"

meta = pd.read_csv(meta_path, delimiter='\t')
meta = meta.rename(columns = {'#SampleID':'SampleID'})
nrows = meta.shape[0]

env_dict = dict()

for i in range(nrows):

    sample = meta.SampleID[i]
    env = meta.empo_1[i]
    env_dict[sample] = env

# print(env_dict)
with open(sbs_path) as r:
    line1 = r.readline()
    line2 = r.readline().strip()
    # line3 = r.readline().strip()
    otus = r.readlines()

# print(line2.split("\t")[0:20])
sites = line2.split("\t")[1:-1] # dont include tax information

otu_table = dict()

site_labels = np.array([ -1 if env_dict[site]=="Free-living" else 1 for site in sites ]) # -1 if freeliving, 1 if host associated

# test
# parts = line3.split("\t")
# otu = parts[1:-1]
# otu_v = np.array(otu, dtype = 'f')

# env_type = site_labels * np.array(otu, dtype='f')
# print((env_type > 0).all())

part_proc = partial(map_proc, site_labels=site_labels)

pool = mp.Pool(mp.cpu_count())

print("%s Processors\n" % mp.cpu_count() )

envs = pool.map(part_proc, otus, chunksize=1)

pool.close()

df = pd.DataFrame(envs, columns = ['ha', 'majority_ha', 'obligate_ha', 'fa', 'majority_fa', 'obligate_fa'])
n = df.shape[0]

ha = np.sum(df.ha)/n
m_ha = np.sum(df.majority_ha)/n
o_ha = np.sum(df.obligate_ha)/n
fa = np.sum(df.fa)/n
m_fa = np.sum(df.majority_fa)/n
o_fa = np.sum(df.obligate_fa)/n

print('host-associated:\t%s' % ha)
print('majority ha:\t%s' % m_ha)
print('obligate ha:\t%s' % o_ha)
print()
print('free-living:\t%s' % fa)
print('majority fa:\t%s' % m_fa)
print('obligate fa:\t%s' % o_fa)