import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random
from scipy.interpolate import interp1d
import time

global JT, hT, Lx, Ly, Zb, Nb, Ns, Iter, config, inc, values
JT0 = 0.5 #J/T ratio
hT = 0.0 #H/T ratio
Lx = 20 #size of lattice in x and y directions
Ly = 20

Zb = 1000 #number of iterations in one block
Nb = 10 #number of blocks
Ns = 1 #number of superblocks

Iter = Zb*Nb*Ns

def initialize():
    config = np.ones((Lx, Ly))
    inc = np.zeros((Lx, Ly), dtype = bool)
    values = np.zeros(Iter)
    return [config, inc, values]

ic = initialize()

config = ic[0] #spins configuration
inc = ic[1] #array used for Wolff cluster growth algorithm, shows whether a spin is included in the cluster
values = ic[2] #array to keep track of total M values, used for error and autocorrelation estimates

def Metropolis_Step():
    x = sp.random.randint(0, Lx-1)
    y = sp.random.randint(0, Ly-1)
    dE = 0.0
    neighbor = x - 1
    if (neighbor > -1):
        dE+=config[neighbor][y]
    neighbor = x + 1
    if (neighbor < Lx):
        dE+=config[neighbor][y]
    neighbor = y - 1
    if (neighbor > -1):
        dE+=config[x][neighbor]
    neighbor = y + 1
    if (neighbor < Ly):
        dE+=config[x][neighbor]

    dE = 2.0*(dE*JT + hT)*config[x, y]
    R = np.exp(-dE)
    if R >= 1.0 or random.random() < R:
        config[x][y] = -config[x][y]

#Metropolis algorithm
def Metropolis():
    A = 0.0
    i = 0
    while i < Iter:
        Metropolis_Step()
        values[i] = abs(sum(sum(config)))/(Lx*Ly)
        A+=abs(sum(sum(config)))
        i+=1
    return A/Iter/(Lx*Ly)

#Wolff algorithm
def Wolff(JT):
    A = 0.0
    i = 0
    random.seed()
    while i<Iter:
        update_Cluster(JT)
        values[i] = abs(sum(sum(config)))/(Lx*Ly)
        A+=abs(sum(sum(config)))
        i+=1
    return A/Iter/(Lx*Ly)


def update_Cluster(JT):
    x = sp.random.randint(0, Lx-1)
    y = sp.random.randint(0, Ly-1)
    st = config[x][y]
    grow_Cluster(JT, st, x, y)
    inc.fill(False)


#Wolff cluster growth algorithm, assuming periodic boundary conditions
def grow_Cluster(JT, st, xi, yi):
    inc[xi][yi] = True
    config[xi][yi] = -config[xi][yi]

    if inc[(xi+1) % Lx][yi]==False and config[(xi+1) % Lx][yi]==st and np.exp(-2.0*JT)<sp.random.random():
        grow_Cluster(JT, st, (xi+1) % Lx, yi)
    if inc[(xi-1) % Lx][yi]==False and config[(xi-1) % Lx][yi]==st and np.exp(-2.0*JT)<sp.random.random():
        grow_Cluster(JT, st, (xi-1) % Lx, yi)
    if inc[xi][(yi+1) % Ly]==False and config[xi][(yi + 1) % Ly]==st and np.exp(-2.0*JT)<sp.random.random():
        grow_Cluster(JT, st, xi, (yi+1) % Ly)
    if inc[xi][(yi-1) % Ly]==False and config[xi][(yi - 1) % Ly]==st and np.exp(-2.0*JT)<sp.random.random():
        grow_Cluster(JT, st, xi, (yi-1) % Ly)



def Error():
    avg = sum(values)/len(values)
    error = 0.0
    for i in range(Ns):
        error += (sum(values[i*Nb*Zb:(i+1)*Nb*Zb])/(Nb*Zb) - avg)**2
    error = np.sqrt(error)/Ns
    return error

def Autocorrelation():
    blockAv = np.zeros(Nb*Ns)
    acf = np.zeros(Nb*Ns)
    for i in range(Nb*Ns):
        blockAv[i] = sum(values[i*Zb:(i+1)*Zb])
    avg = sum(blockAv)/len(blockAv)

    for j in range(Nb*Ns):
        num = sum((blockAv[0:(Nb*Ns - j)] - avg)*(blockAv[j:Nb*Ns] - avg))
        den = sum((blockAv[0:(Nb*Ns - j)] - avg)**2)
        acf[j] = num/den

    return acf


#used for plotting a magnetization curve
def curve(JTs):
    vecw = np.vectorize(Wolff)
    y = vecw(JTs)
    return y
