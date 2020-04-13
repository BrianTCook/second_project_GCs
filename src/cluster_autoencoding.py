#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 16:45:15 2020

@author: BrianTCook
"""

import pandas as pd
import numpy as np
from keras.layers import Input, Dense
from keras.models import Model
from keras import optimizers
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

def get_autoencoded_result(Norbiters, ndim_input):

    #load and prepare data
    datadir = '/Users/BrianTCook/Desktop/Thesis/second_project_GCs/data/'
    cluster_populations = np.loadtxt(datadir + 'Nstars_in_clusters.txt')
    cluster_populations = list(cluster_populations[:Norbiters])
    
    labels = []
    
    for logN in range(int(np.log2(Norbiters))+1):
        labels += [ logN for i in range(int(cluster_populations[logN])) ]
    
    data_filename = glob.glob(datadir+'enbid_files/*_00400_Norbiters_%i.ascii'%(Norbiters))
    
    data_6D = np.loadtxt(data_filename[0])
    data_3D = data_6D[:, 0:3]
    data_2D = data_6D[:, 1:3]    
    
    vectors = input_vecs[:,2:]
    labels = input_vecs[:,:2]
    
    x_train, x_test, y_train, y_test = train_test_split(vectors, labels, test_size=0.33, random_state=1)
    
    #Setup autoencoder structure
    input_vector = Input(shape=(16,)) #Use 16 floats vector as input
    encoded_1 = Dense(500, activation='selu')(input_vector)
    encoded_2 = Dense(250, activation='selu')(encoded_1)
    middle_layer = Dense(2, activation='elu')(encoded_2)
    decoded_1 = Dense(250, activation='selu')(middle_layer)
    decoded_2 = Dense(500, activation='selu')(decoded_1)
    decoded_output = Dense(16, activation='sigmoid')(decoded_2)
    
    '''
        Note: I have been playing around with all kinds of activation
        functions above, as well as optimizers and loss functions below.
        Please share your thoughts on which would be reasonable to use.
        
        Nadam and rmsprop seem to have the most promising results.
    '''
    
    autoencoder = Model(input_vector, decoded_output)
    autoencoder.compile(optimizer='rmsprop', loss='mae')
    autoencoder.summary()
    
    #Train the autoencoder
    autoencoder.fit(x_train, x_train,
                    epochs=30, #also used many different epoch values/batch sizes
                    batch_size=4,
                    shuffle=True,
                    validation_data=(x_test, x_test))
    
    autofitted = autoencoder.predict(vectors)
    
    #compare to estimate performance
    autofitted[1]
    x_test[1]
    
    #separate encoder
    encoder = Model(input_vector, middle_layer)
    encoder.compile(optimizer='rmsprop', loss='mae')
    
    #predict encoded vectors from input vectors
    encoded_vecs = encoder.predict(vectors)
    
    #separate decoder
    encoded_vector = Input(shape=(2,))
    deco1 = autoencoder.layers[-3](encoded_vector)
    deco2 = autoencoder.layers[-2](deco1)
    decoout = autoencoder.layers[-1](deco2)
    decoder = Model(encoded_vector, decoout)
    
    #predict decoded vectors from encoded vectors
    decoded_vecs = decoder.predict(encoded_vecs)
    
    #Confirm that separate encoding/decoding does the same as the autoencoder
    np.array_equal(autofitted, decoded_vecs)
    
    #Write labels with corresponding encoded values to file
    output = np.concatenate((labels,encoded_vecs), axis=1)
    np.savetxt("encoded_with_labels.txt", output, fmt='%s', delimiter=',')
    
    return 0