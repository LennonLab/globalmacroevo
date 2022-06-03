# MicroSpeciation

Macroevolutionary modeling of bacterial and archaeal evolution. 

16S rRNA substitution rates taken from [Kuo and Ochman (2009)](https://biologydirect.biomedcentral.com/articles/10.1186/1745-6150-4-35).

### Figure code

Figure 1: `/code/bd_expectations.R`
Figure 2: `/code/mass_extinction.R`
Figure 1S: `/code/probability_per_magnitude.R`
Figure 2S: `/code/CladeRatesMassExtinction.R`

### EMP analysis

Code that takes Earth Microbiome Project (EMP) data on species and sites to determine global proportions of host-associated and free-living taxa. 

EMP metadata: `/data/EMP/emp_qiime_mapping_release1_20170912.tsv`
EMP biom table:  https://zenodo.org/record/890000
EMP job script: `/code/emp.script`
EMP analysis: `/code/processEMB.py`
Output: `/data/EMP/emp_ratio.txt`

### Other scripts

`/code/functions.R`
Data transformation and utility functions for use in other scripts.

`/code/speciation_from_sub.R`
Calculation of speciation rates from Kuo and Ochman (2009) 16S rRNA substitution rates. 

## ancient-seq

Attempting to use ancient sequences and present-day sequences to directly estimate 16S rRNA substitution rates. 

Contains ancient sequences used and code used to make phylogenies.

Ancient sequences appear to highly damaged and unable to be used for substitution rate calculation. 