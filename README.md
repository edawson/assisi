assisi: another somatic signature simulator
-------------------------------------------
Eric T. Dawson  
September 27 2017


## What is it?
`assisi` is a python package that simulates mutational signature profiles (i.e. counts / proportions of trinucleotide context mutations) using Monte Carlo simulation.

## Required / useful files
We've included the cosmic mutational signatures probabilities from deconstructSigs (https://github.com/raerose01/deconstructSigs).

## Running assisi
`python assisi.py -h` for help.  

## Typical usage: `
The following command will generate a single mutational catalog with 1000 total mutations.
`-r` specifies that assisi should generate a random number of signatures (1-10) and assign
random proportions to each of those ten signatures. The sum of the proportions is 1.0,
with some error.  
`-i` is the COSMIC SBS96 mutational context probabilities; each row is a signature, and each column
is the proportion of that mutational context in the signature.

```
python assisi.py -i cosmic_file.tsv -n 1000 -r > catalog.tsv
```

Rules for COSMIC sig files:  
    No headers  
    Tab separated  
    No restriction on number of rows  
    There must be 96 columns (for now)  


The next version of assisi will use a tidy format for this data.

The output of assisi when run as above is a single line with the following columns:  
```
<number of signatures> <Python dictionary of signature : proportion> <first context count> ... <last context count>
```

Multiple catalogs can be generated by passing the `-m` flag and a desired number of catalogs:
```
## Simulate 1000 catalogs
python assisi.py -i cosmic_file.tsv -n 1000 -r -m 1000 > catalogs.tsv
```

Assisi uses the Python multiprocessing package to optionally generate catalogs in parallel.
We recommend using this feature when simulating more than one catalog, as it decreases runtime
proportional to the number of available compute cores and does not affect results in any way.
```
## Simulate 1000 catalogs, with 4 threads
python assisi.py -i cosmic_file.tsv -n 1000 -r -m 1000 -t 4 > catalogs.tsv
```


## Installation
`assisi` is written in barebones python3 and is backwards-compatible with python 2.7 - no installation necessary.

