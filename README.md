assisi: another somatic signature simulator
-------------------------------------------
Eric T. Dawson  
September 27 2017


## What is it?
`assisi` is a python package that simulates mutational signature profiles (i.e. counts / proportions of trinucleotide context mutations) using Monte Carlo simulation.

## Required / useful files
We've included the cosmic mutational signatures probabilities from deconstructSigs (https://github.com/raerose01/deconstructSigs).

## Running assisi
`python assisi.py /path/to/mutational/props/file.tsv`  

Rules:  
    No headers  
    Tab separated  
    No restriction on number of rows  
    There must be 96 columns (for now)  

## Installation
    `assisi` is written in barebones python2.7 - no installation necessary.

