import os
import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt

#================
# Simulation parameters
t = 0.0
tmax = 1
dt = 0.01
tstepsNb = tmax/dt
YM = 1


#================
# Initial parameters
l0i = np.zeros(shape=(int(tstepsNb+1),4))
l0i[0,:] =  [10,11,12,13]
yi = [-1.00, -0.33, 0.33, 1.00]


#================
# Output Variables
li = np.zeros(shape=(int(tstepsNb+1),4))
li[0, :] = [np.mean(l0i[0,:])]*len(l0i[0,:])
# l = np.zeros(shape=(int(tstepsNb+1)))
ki = np.zeros(shape=(int(tstepsNb+1),4))
# k = np.zeros(shape=(int(tstepsNb+1)))


def equations(params):
    liPrev = params[0:4]
    # print("liPrev", liPrev)
    l0i = params[4:8]
    # print("l0i", l0i)
    kiPrev = params[8:12]
    # print("ki", kiPrev)
    yi = params[12:16]
    # print("yi", yi)
    # , l0i, kiPrev, yi = params
    liNext = sum( ( liPrev[i] * ( 1 + ki[i] * yi[i] ) - l0i[i] ) * (liPrev[i] * yi[i]) for i in range(0, len(yi)) )
    # print(len(liNext))
    kiNext = sum( ( liPrev[i] * ( 1 + ki[i] * yi[i] ) - l0i[i] ) * ( 1 + ki[i] * yi[i] ) for i in range(0, len(yi)) )
    # print(len(kiNext))
    sol = np.concatenate((liNext, l0i), axis = 0)
    sol = np.concatenate((sol , kiPrev), axis = 0)
    sol = np.concatenate((sol , yi), axis = 0)
    # print(sol)
    return sol

def stress(li, l0i, YM):
    stress = YM * (li - l0i)
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
    tVec.append(t)
    tIndex += 1

    params = [li[tIndex-1,:], l0i[tIndex-1,:], ki[tIndex-1,:], yi]
    # print(params)
    sol = fsolve(equations, params )
    # print(sol)
    li[tIndex, :] = sol[0:4]
    ki[tIndex, :] = sol[8:12]

    st = stress(li[tIndex-1,:], l0i[tIndex-1,:], 1)
    # print(st)

    gr = growth(st, l0i[tIndex-1], dt)
    # print(gr)

    l0i[tIndex] = gr

# print(li)
print(li)
plt.figure()
plt.plot(tVec, ki[:,0])
plt.show()
