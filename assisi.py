import sys
import random
import multiprocessing as mp
import argparse


def print_err(s):
    sys.stderr.write(s)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type = str,
     help = "Input mutational proportions / counts file.", dest = "infile")
    parser.add_argument("-d", dest="sigdict", type = str,
    help = "A signature:proportion dictionary file, CSV format, one [sig,prop] pair per line")
    parser.add_argument("-r", dest="random", action = "store_true", default = False,
        help = "Randomly initialize sigs / amounts. Default: False")
    parser.add_argument("-c", action="store_true",
     help = "Input is raw counts, rather than proportions. Default: False", dest="iscounts", default = False)
    parser.add_argument("-n", type = int, dest="num",
     help = "The number of mutations to simulate. Default: 1000", default = 1000)
    parser.add_argument("-m", type = int,
     dest="numsamples", help = "The number of mutational spectra to simulate. Default: 1", default = 1)
    parser.add_argument("-t", type = int, dest="threads", default = 1,
     help = "The number of multiprocessing threads to use, parallelized by spectrum. Default: 1")
    parser.add_argument("-e", dest="eps", type = float,
     help = "The proportion of error to incorporate. Default 0.0", default = 0.0)
    parser.add_argument("-p", action = "store_true", dest = "showdist",
     default = False, help = "Print distribution as ugly ASCII to terminal.")

    return parser.parse_args()

def parse_sig_dict(dict_file):
    d = {}
    with open(dict_file, "r") as ifi:
        for line in ifi:
            tokens = line.strip().split(",")
            if len(tokens) > 3:
                print "ERROR: too many columsn in sig_dict file."
                print "There should be only two values, separated by a comma, per line."
                print "E.g. \"1,0.5\""
                exit(9)
            d[int(tokens[0].strip())] = float(tokens[1].strip())
    return d


def parse_sigs(sig_file):
    siggy = []
    with open(sig_file, "r") as s:
        for line in s:
            tokens = line.strip().split("\t")
            tokens = [float(i) for i in tokens]
            siggy.append(tokens)
    return siggy

def rmse(sig, beta):
    err = 0.0
    for i in len(sig):
        err += (sig[i] - beta[i])^2
    return err

## Args: a 96-length vector of mutation counts
def counts_to_props(counts):
    total = sum(counts)
    p = [(float(i) / float(total)) for i in counts]
    return p

## Args:
##  probs: A 96-length vector of probabilities for that sample
##  index: the index of the desired mut category
##  value: the value from the random number generator
def sample_prob(probs, index, val):
    if val <= probs[index]:
        return True
    return False

## Args:
##  sig_vec: an empty vector to fill with mutation counts
def sample_sig(probs, number):
    ## The following line caches a function in python,
    ## and since we call this one many times, we want to do that.
    csp = sample_prob
    ret = [0 for i in xrange(0, len(probs))]
    count = 0
    last = 0
    while count < number:
        x = random.randint(0, len(probs) -  1)
        val = random.random()
        if csp(probs, x, val):
            ret[x] += 1
            last = x
            count += 1
    return ret, last 

## Args:
##  sig_and_amt: a dictionary of the signature index
##      and the amount of that signature as a percentage of total mutations ( a float.)
##  number: the total number of mutations to produce
##  probs: an array (list of lists), numSigs x numContexts, containing the probabilities for each
##  mutation type within each signature
def sample_sig_list(sig_and_amt, number, probs, num_sigs = 30, eps = 0.0):
    if len(sig_and_amt) == 0 and eps < 0.01:
        print "FAILING: low epsilon and no signature proportions specified"
        print "Please specify either of these values to prevent running infinitely."
        exit(9)

    cspp = sample_prob
    sigs = [i for i in xrange(0, num_sigs)]
    amts = [0.0 for i in xrange(0, num_sigs)]
    for i in sig_and_amt:
        amts[i - 1] = float(sig_and_amt[i])
    
    ## Amounts is a now a probability vector
    ## to introduce process noise, we'll make the minimum probability
    ## of any signature vary by eps
    amts = [eps if i == 0 else i for i in amts]

    ## Now we're going to sample the signatures
    ret = [0.0 for i in xrange(0, 96)]
    count = 0
    while count < number:
        x = random.randint(0, len(sigs) - 1 )
        val = random.random()
        if cspp(amts, x, val):
            dist, last = sample_sig(probs[x], 1)
            ret[last] += 1
            count += 1
    return ret

def print_dist(sig, use_tight_ranges = False):
    ostr = ""
    ranges = []
    if use_tight_ranges:
        ranges = [0.99,0.98,0.97,0.96,0.95,0.94,
        0.93,0.92,0.91,0.90,0.89,0.88,0.87,0.86,
        0.85,0.84,0.83,0.82,0.81,0.80,0.79,0.78,
        0.77,0.76,0.75,0.74,0.73,0.72,0.71,0.70,
        0.69,0.68,0.67,0.66,0.65,0.64,0.63,0.62,
        0.61,0.60,0.59,0.58,0.57,0.56,0.55,0.54,
        0.53,0.52,0.51,0.50,0.49,0.48,0.47,0.46,
        0.45,0.44,0.43,0.42,0.41,0.40,0.39,0.38,
        0.37,0.36,0.35,0.34,0.33,0.32,0.31,0.30,
        0.29,0.28,0.27,0.26,0.25,0.24,0.23,0.22,
        0.21,0.20,0.19,0.18,0.17,0.16,0.15,0.14,
        0.13,0.12,0.11,0.10,0.09,0.08,.07,0.6,0.05,
        0.04,0.03,0.02,0.01, 0.005]
    else:
        ranges = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6,
        0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
    
    m = max(sig)
    for i in ranges:
        if i > m:
            continue
        ostr += '{0:05.1%}'.format(i) + " "
        for s in sig:
            if s >= i:
                ostr += "|"
            else:
                ostr += " "
            ostr += " "
        ostr += "\n"
    ostr += "======"
    for s in sig:
        ostr += "-="
    print_err(ostr)
    return ostr


if __name__ == "__main__":
    #test = [random.random() for i in xrange(0, 96)]
    #print_dist(test)
    #print test
    args = parse_args()

    probs = parse_sigs(args.infile)
    if args.iscounts:
        probs = [counts_to_props(i) for i in probs]


    n_muts = args.num
    n_samples = args.numsamples

    if n_samples > 1 or args.threads > 1:
        print "threading not implemented"
        exit()
        p = mp.pool(args.threads)


    sigm = {}

    if args.sigdict is not None:
        sigm = parse_sig_dict(args.sigdict)

    elif args.random:
        rem = 1.0
        while rem > 0.01:
            i = random.randint(0, 30)
            if i not in sigm:
                nex = random.uniform(0, rem)
                sigm[random.randint(0,30)] = nex
                rem = rem - nex
        print "Randomly simulating with", len(sigm), "signatures"
        print sigm
    else:
        ## Enforce flat probabilities
        for i in xrange(0, len(probs)):
            sigm[i] = 1.0

    
    d = sample_sig_list(sigm, n_muts, probs, len(probs), args.eps)
    c = counts_to_props(d)
    
    if (args.showdist):
        print_dist(c, True)
    print "\t".join([str(int(i)) for i in d])