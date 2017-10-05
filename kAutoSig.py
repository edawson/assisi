import argparse
import sys
import math
import numpy as np
from keras.layers import Input, Dense
from keras.models import Model
from assisi import print_dist, counts_to_props


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest = "infile", type = str, required = True, help = "Input sample X mutation counts matrix")
    return parser.parse_args()

def prep_data(dat):
    ret = None
    return ret

def train(sigs, enc_dim):
    inpt = Input(shape=(96,))
    encoded = Dense(enc_dim, activation='relu')(inpt)
    decoded = Dense(96, activation='sigmoid')(encoded)

    autoencoder = Model(inpt, decoded)

    encoder = Model(inpt, encoded)

    enc_inp = Input(shape = (enc_dim,))

    dec_layer = autoencoder.layers[-1]

    decoder = Model(enc_inp, dec_layer(enc_inp))

    autoencoder.compile(optimizer='adadelta', loss='binary_crossentropy')



    autoencoder.fit(xtr, xtr,
                epochs=50,
                batch_size=25,
                shuffle=True,
                validation_data=(xtr, xtr))

    en_d = encoder.predict(xtr)
    de_d = decoder.predict(en_d)

    print xtr[1]
    print_dist(xtr[1], True)
    print de_d[1]
    print_dist(de_d[1], True)
    return

def test():
    return

if __name__ == "__main__":

    args = parse_args()

    sigs = []
    with open(args.infile, "r") as ifi:
        for line in ifi:
            sigs.append(counts_to_props([int(i) for i in line.strip().split("\t")]))
    
    #xtr = [np.array(i) for i in sigs]
    xtr = np.array(sigs)

    enc_dim = 30
    
    train(sigs, enc_dim)


