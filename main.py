import os
import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt

#================
# Simulation parameters
t = 0.0
tmax = 0.1
dt = 0.01
tstepsNb = tmax/dt
YM = 1.0


#================
# Initial parameters
l0i = np.zeros(shape=(int(tstepsNb+2),4))
l0i[0,:] =  [10,11,12,13]
yi = [-1.00, -0.33, 0.33, 1.00]


#================
# Output Variables
l = np.zeros(shape=(int(tstepsNb+2)))
l[0] = np.mean(l0i[0,:])
# l = np.zeros(shape=(int(tstepsNb+1)))
k = np.zeros(shape=(int(tstepsNb+2)))
k[0] = 1


def equations(params):
    lPrev, kPrev = params
    lNext = sum( ( lPrev * ( 1 + k * yi[i] ) - Fl0i[i] ) * (lPrev * yi[i]) for i in range(0, len(yi)) )
    kNext = sum( ( lPrev * ( 1 + k * yi[i] ) - Fl0i[i] ) * ( 1 + kPrev * yi[i] ) for i in range(0, len(yi)) )
    print("shape", kNext.shape)
    return lNext[tIndex-1], kNext[tIndex-1]

def stress(l, l0i, YM):
    stress = YM * (l - l0i)
    return stress

def growth(stress, l0i, dt):
    stressMinus = [float(-1*s) for s in stress]
    # print(stress)
    func = 1/(1+np.exp(stressMinus))
    Grow = []
    for i in range(0, len(stress)):
        # print(func[i], l0i[i])
        Grow.append(l0i[i] + (func[i] * l0i[i])*dt)
    return Grow


tIndex = 0
tVec = [t]
while t < tmax:
    t += dt
    print("t = ", t)
    tVec.append(t)
    tIndex += 1

    params = [l[tIndex-1], k[tIndex-1]] #  [li[tIndex-1,:], l0i[tIndex-1,:], ki[tIndex-1,:], yi]
    # print(yi)
    Fl0i = l0i[tIndex-1,:]
    print(Fl0i)
    # globParams = np.concatenate((l0i[tIndex-1],yi), axis=0)
    # print(globParams)
    # print(params)
    l[tIndex], k[tIndex] = fsolve(equations, params )
    # print(sol)
    # l[tIndex, :] = sol[0:4]
    # k[tIndex, :] = sol[8:12]
    #
    st = stress(l[tIndex-1], l0i[tIndex-1,:], YM)
    print("stress", st)
    #
    gr = growth(st, l0i[tIndex-1], dt)
    print("growth", gr)
    #
    l0i[tIndex] = gr

# print(li)
print(l)
plt.figure()
plt.plot(tVec, l[:])
plt.show()
