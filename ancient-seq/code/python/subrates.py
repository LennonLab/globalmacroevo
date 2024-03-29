## subrates.py
## Ford Fishman
## 7/6/21

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-wd',  type=str, help="Path to directory with input files")
parser.set_defaults(single=True)
arguments = parser.parse_args()
wd = arguments.wd

age_file = wd + "accessions_age.csv"
div_file = wd + "divergence_silva.csv"
output = wd + "subrates_1My.csv"

ages = dict()

with open(age_file, "r") as f:

    f.readline()

    for line in f:



        line = line.strip()
        parts = line.split(',')

        acc = str(parts[0])
        # print(acc)

        if not acc[-2]=='.': # some accessions are missing trailing decimal
            acc += '.1'
        
        ages[acc] = {'low':float(parts[1]), 'high': float(parts[2])}
        
sub_rates = dict()

with open(div_file, "r") as r:
    
    with open(output, 'w') as f:

        f.write('%s,%s,%s,%s,%s,%s\n' % ('acc', 'low_sub','high_sub','age_low', 'age_high', 'div'))
    

        for line in r:

            line = line.strip()
            parts = line.split(',')

            acc = str(parts[0])
            div = float(parts[1])

            age_low = ages[acc]['low'] 
            age_high = ages[acc]['high'] 

            if age_high < 1:
                continue

            if age_low == 0:
                high_sub = 0
            
            else: 
                high_sub = div/age_low
                
            low_sub = div/age_high

            sub_rates[acc] = {'low': low_sub, 'high': high_sub}

            if (high_sub == 0 and low_sub == 0) or (np.isnan(high_sub) and np.isnan(low_sub)):
                continue

            f.write('%s,%s,%s,%s,%s,%s\n' % (acc, low_sub, high_sub, age_low, age_high, div))




