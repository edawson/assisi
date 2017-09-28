import sys
from assisi import rmse
from assisi import counts_to_props 
import math

def recount(sig, num):
    return [math.floor(i * num) for i in sig]

def diff(sig, beta):
    s = 0.0
    for i in xrange(0, len(sig)):
        s += abs(sig[i] - beta[i])
    return s

if __name__ == "__main__":

    num = 10000
    sigs = []
    rands = []
    with open(sys.argv[1], "r") as ifi:
        for line in ifi:
            tokens = line.strip().split("\t")
            tokens = [int(i) for i in tokens]
            tokens = counts_to_props(tokens)
            sigs.append(tokens)
    
    with open(sys.argv[2], "r") as ifi:
        for line in ifi:
            tokens = line.strip().split("\t")
            tokens = [int(i) for i in tokens]
            tokens = counts_to_props(tokens)
            rands.append(tokens)

    for i in xrange(0, len(sigs)):
        r = rmse(sigs[i], rands[i])
        d = diff(sigs[i], rands[i])
        print r, d
