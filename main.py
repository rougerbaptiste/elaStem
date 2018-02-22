import os
import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt

#================
# Simulation parameters
t = 0.0
tmax = 10.0
dt = 0.01
tstepsNb = tmax/dt
YM = 1.0
fiberNb = 4
prodRate = 0.8
degRate = 2.0
hormC = [0.5] * fiberNb
print(hormC)

#================
# Initial parameters
l0i = np.zeros(shape=(int(tstepsNb+1),fiberNb))
l0i[0,:] =  np.linspace(10,10+fiberNb, fiberNb)
yi = np.linspace(-1, 1, fiberNb) # [-1.00, -0.33, 0.33, 1.00]



#================
# Output Variables
l = np.zeros(shape=(int(tstepsNb+1)))
l[0] = np.mean(l0i[0,:])
# l = np.zeros(shape=(int(tstepsNb+1)))
k = np.zeros(shape=(int(tstepsNb+1)))
k[0] = 1


def equations(params):
    lPrev, kPrev = params
    lNext = sum( ( lPrev * ( 1 + kPrev * yi[i] ) - Fl0i[i] ) * (lPrev * yi[i]) for i in range(0, len(yi)) )
    kNext = sum( ( lPrev * ( 1 + kPrev * yi[i] ) - Fl0i[i] ) * ( 1 + kPrev * yi[i] ) for i in range(0, len(yi)) )
    return lNext, kNext

def stress(l, l0i, YM):
    stress = YM * (l - l0i)
    return stress

def growth(stress, l0i, dt):
    stressMinus = [float(-1*s) for s in stress]
    # print(stress)
    func = 1/(1+np.exp(stressMinus))
    Grow = []
    for i in range(0, len(stress)):
        hormC[i] = hormC[i] * (1 - dt/degRate) + max(dt * func[i], 0)
        # hormC[i] = hormC[i] * (1 - dt/degRate) + max(dt * prodRate * stress[i], 0)
        Grow.append(l0i[i]* (1 + max(dt * hormC[i],0)))
    print(hormC)
    print(Grow)
    return Grow


tIndex = 0
tVec = [t]
while t < tmax:
    t += dt
    if(t > tmax): break
    print("t = ", t)
    tVec.append(t)
    tIndex += 1

    params = [l[tIndex-1], k[tIndex-1]] #  [li[tIndex-1,:], l0i[tIndex-1,:], ki[tIndex-1,:], yi]
    # print(yi)
    Fl0i = l0i[tIndex-1,:]
    # globParams = np.concatenate((l0i[tIndex-1],yi), axis=0)
    # print(globParams)
    # print(params)
    l[tIndex], k[tIndex] = fsolve(equations, params )
    # print(sol)
    # l[tIndex, :] = sol[0:4]
    # k[tIndex, :] = sol[8:12]
    #
    st = stress(l[tIndex-1], l0i[tIndex-1,:], YM)
    #
    gr = growth(st, l0i[tIndex-1], dt)
    # print(gr)
    #
    l0i[tIndex, :] = gr

# print(l)
#
# print(k)
# print(len(tVec), k.shape)
plt.figure()
plt.plot(tVec, k[:])
plt.xlabel('time')
plt.ylabel('curvature')
plt.figure()
plt.plot(tVec, l)
plt.xlabel('time')
plt.ylabel('length')
# plt.axis((0,2,10,14))
plt.show()
