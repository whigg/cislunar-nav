'''
TODO:
 - get rid of user position, all we care about is change in position and velocity
 - if that doesn't work, switch to evaluating force
'''

from Parsing import *

import os
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor

def parse_and_scale():
    dir = 'data/training/'
    files = []
    for file in os.listdir(dir):
        f = os.path.join(dir, file)

        if os.path.isfile(f):
            files.append(f)
    
    for i, f in enumerate(files):
        if i == 0:
            pos, vel, moon, earth, sun, force, time = parseGmatData(f, gmatReport=True)
        else:
            pt, vt, mt, et, st, ft, tt = parseGmatData(f, gmatReport=True)
            pos._vec = np.concatenate((pos._vec, pt._vec), axis=0)
            vel._vec = np.concatenate((vel._vec, vt._vec), axis=0)
            moon._vec = np.concatenate((moon._vec, mt._vec), axis=0)
            earth._vec = np.concatenate((earth._vec, et._vec), axis=0)
            sun._vec = np.concatenate((sun._vec, st._vec), axis=0)
            force._vec = np.concatenate((force._vec, ft._vec), axis=0)
            time = np.concatenate((time, tt), axis=0)

    dMoon  = (moon._vec - pos._vec)[:-1,:]
    dEarth = (earth._vec - pos._vec)[:-1,:]
    dSun   = (sun._vec - pos._vec)[:-1,:]
    dt = time[1:] - time[:-1]

    # X = np.concatenate((pos._vec[:-1,:], vel._vec[:-1,:], dMoon, dEarth, dSun, dt), axis=1)
    # Y = np.concatenate((pos._vec[1:,:], vel._vec[1:,:]), axis=1)
    X = np.concatenate((pos._vec, time), axis=1)
    Y = force._vec

    # Scale the data (X only)
    xscaler = StandardScaler()
    xscaler.fit(X)
    X = xscaler.transform(X)

    yscaler = StandardScaler()
    yscaler.fit(Y)
    Y = yscaler.transform(Y)

    return X, Y, xscaler, yscaler

def train_model(X, Y):
    Neural = MLPRegressor(hidden_layer_sizes=(150,50), solver='adam', random_state=69, max_iter=1500,
                          tol=0.00001, learning_rate='adaptive', verbose=True).fit(X, Y)

    return Neural

