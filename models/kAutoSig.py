import argparse
import sys
import math
import numpy as np
from keras.layers import Input, Dense
from keras.models import Model
from keras import backend as K
from assisi import print_dist, counts_to_props


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest = "infile", type = str, required = True, help = "Input sample X mutation counts matrix")
    return parser.parse_args()

def prep_data(dat):
    ret = np.array(dat)
    return ret

def denoise(dat):
    return

def train(sigs, enc_dim):
    inpt = Input(shape=(96,))
    encoded = Dense(enc_dim, activation='relu')(inpt)
    #encoded = Dense(60, activation='relu')(inpt)
    #encoded = Dense(30, activation='relu')(encoded)
    #encoded = Dense(enc_dim, activation='relu')(encoded)
    decoded = Dense(96, activation='sigmoid')(encoded)

    autoencoder = Model(inpt, decoded)

    encoder = Model(inpt, encoded)

    enc_inp = Input(shape = (enc_dim,))

    dec_layer = autoencoder.layers[-1]

    decoder = Model(enc_inp, dec_layer(enc_inp))

    autoencoder.compile(optimizer='adadelta', loss='binary_crossentropy')

    grab_layer = K.function([encoder.layers[0].input],
                                  [encoder.layers[0].output])


    autoencoder.fit(xtr, xtr,
                epochs=100,
                batch_size=30,
                shuffle=True,
                validation_data=(xtr, xtr))

    en_d = encoder.predict(xtr)
    de_d = decoder.predict(en_d)

    print xtr[0]
    print_dist(xtr[0], True)
    print de_d[0]
    print_dist(de_d[0], True)

    #print grab_layer([xtr])[0][0]
    # for i in grab_layer([xtr[0]])[0]:
    #     print i
    return

def validate():
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
    xtr = prep_data(sigs)

    enc_dim = 5
    
    train(sigs, enc_dim)


