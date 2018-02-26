import os
import numpy as np
from scipy.optimize import fsolve
import math
import matplotlib.pyplot as plt

#
# (S) means the varaible is a scalar
# (V) means the varaible is a vector
# (M)(x,y) means the varaible is a matrix of size x by y
#

degRate = prodRate / ratio               # (S) Degradation rate of the hormone of growth

ProdDegRatio = [0.2]

outputMatrixL = []
outputMatrixK = []

for ratio in ProdDegRatio:
    print("\n\n=============" + str(ratio) + "=============== \n")
    #================
    # Simulation parameters
    t = float(0.0)                     # (S) Starting time of the simulation
    tmax = float(50.0)                 # (S) Ending time of the simulation
    dt = float(0.01)                   # (S) Time step
    tstepsNb = tmax/dt          # (S) Number of steps of the simulation
    YM = 5.0                    # (S) Young Modulus : elasticity coefficient of the stem
    fiberNb = 4                 # (S) Number of fibers/spring to consider
    prodRate = 1.0              # (S) Production rate of the hormone of growth
    degRate = 5.0               # (S) Degradation rate of the hormone of growth
    hormC = [0.5] * fiberNb     # (V) Concentration in hormone in the fibers
    hormC2 = [0.5] * fiberNb     # (V) Concentration in hormone in the fibers

    #================
    # Initial parameters
    l0i = np.zeros(shape=(int(tstepsNb+1),fiberNb)) # (M)(simulation length, fiber number) "Original" length of the fibers
    l0i[0,:] =  np.linspace(10,10+fiberNb, fiberNb) # initialisation of the above
    yi = np.linspace(-1, 1, fiberNb)                # (V) position of the fibers in the curve


    #================
    # Output Variables
    l = np.zeros(shape=(int(tstepsNb+1)))           # (V) Actual length of the stem
    l[0] = np.mean(l0i[0,:])                        # Initialisation of the above with mean value of l0i
    k = np.zeros(shape=(int(tstepsNb+1)))           # (V) Curvature of the stem
    k[0] = 1                                        # Initialisation of the above


    def equations(params):
        """
        This function contains the equations to solve to get the new length and curvature of the stem.
        It takes params as an input with inside :
            - (S) lPrev : previous length of the stem
            - (S) kPrev : previous curvature of the stem

        It returns :
            - (S) lNext : the new length of the stem
            - (S) kNext : the new curvature of the stem

        It also needs global variables such as :
            - (V) yi   : the spatial position of each fiber
            - (V) Fl0i : the "original" length of each fiber
        """
        lPrev, kPrev = params
        lNext = sum( ( lPrev * ( 1 + kPrev * yi[i] ) - Fl0i[i] ) * (lPrev * yi[i]) for i in range(0, len(yi)) )         # This computes the partial derivative against l
        kNext = sum( ( lPrev * ( 1 + kPrev * yi[i] ) - Fl0i[i] ) * ( 1 + kPrev * yi[i] ) for i in range(0, len(yi)) )   # This computes the partial derivative against kappa
        return lNext, kNext

    def stress(l, l0i, YM):
        """
        This function computes the stress of each fiber.
        It takes as input :
            - (S) l   : the actual length of the stem at this time step
            - (V) l0i : the current "original" length of each fiber
            - (S) YM  : the Young Modulus that represents the elacticity of the fiber
        It returns :
            - (V) stress : the stress of each fiber of the stem
        """
        stress = YM * (l - l0i)
        return stress

    def growth(stress, l0i, dt, function):
        """
        This function computes the new "original" length of each fiber
        It takes as input :
            - (V) stress : the stress of each fiber of the stem
            - (V) l0i    : the previous "original" length of each fiber
            - (S) dt     : the time step
            - (S) function : a number the specifies the function to use for the growth
                - 0 : classical hormone production + degradation
                - 1 : hormone production following a sigmoide function + degradation

        It returns :
            - (V) l0i : the new "original" length of each fiber

        It also needs global variables such as :
            - (S) degRate  : the degradation rate of the hormone
            - (S) predRate : the production rate of the hormone
        """
        Grow = []                                                                           # (V) Initialisation of the new "original" length of the fibers
        if (function == 0):                                                                 # the user chose function 0 (cf above)
            for i in range(0, len(stress)):                                                 # we compute the growth of each fiber
                hormC[i] = hormC[i] * (1 - dt/degRate) + max(dt * prodRate * stress[i], 0)
                Grow.append(l0i[i]* (1 + max(dt * hormC[i],0)))

        if (function == 1):                                                                 # the user chose function 1 (cf above)
            stressMinus = [float(-1*s) for s in stress]                                     # computes each stress * -1
            func = 1/(1+np.exp(stressMinus))                                                # the sigmoide function of the stress
            for i in range(0, len(stress)):                                                 # we compute the growth for each fiber
                hormC[i] = hormC[i] * (1 - dt/degRate) + max(dt * func[i], 0)
                Grow.append(l0i[i]* (1 + max(dt * hormC[i],0)))

        if (function == 2):
            stressMinus = [float(-1*s) for s in stress]                                     # computes each stress * -1
            func = 1/(1+np.exp(stressMinus))                                                # the sigmoide function of the stress
            for i in range(0, len(stress)):
                hormC[i] = hormC[i] * (1 - dt/degRate) + max(dt * func[i], 0)
                hormC2[i] = hormC2[i] * (1 - dt/degRate) + max(dt * hormC[i], 0)
                Grow.append(l0i[i]* (1 + max(dt * hormC2[i],0)))
        return Grow


    tIndex = 0                                              # (S) the index of the current time step
    tVec = [t]                                              # (V) initialisation of the vector containing the time of each time step
    while t <= tmax:
        t += dt
        if(t > (tmax + dt/100)): break                      # this function is here to avoid having an extra time step

        print("t = " +  str(t))                             # we update the time step at which we are
        tVec.append(t)
        tIndex += 1

        params = [l[tIndex-1], k[tIndex-1]]                 # (V) we put the previous l and k in a vector for the function
        Fl0i = l0i[tIndex-1,:]                              # (V) we put the previous "original" length of the fibers in a vector to sned it to the function
        l[tIndex], k[tIndex] = fsolve(equations, params )   # (V), (V) we get the newly computed values of l and k and store it in the corresponding vectors

        fibersStress = stress(l[tIndex-1], l0i[tIndex-1,:], YM) # (V) we compute the stress of each fibers

        l0i[tIndex, :] = growth(fibersStress, l0i[tIndex-1], dt, 1) # (M) we compute the new "original" lengths of each fiber and store it in the corresponding matrix


    outputMatrixL.append(l)
    outputMatrixK.append(k)

    plt.figure()

    plt.plot(tVec, l)
    plt.legend(ProdDegRatio, loc='upper left')
    plt.suptitle("Length of the stem")
    # plt.savefig(filenameL, dpi=900)

    plt.figure()
    plt.plot(tVec, k)
    plt.legend(ProdDegRatio, loc='upper left')
    plt.suptitle("Curvature of the stem")
    # plt.savefig(filenameK, dpi=900)

    plt.show()

filenameL = '-'.join(str(x) for x in ProdDegRatio) + "L.pdf"
filenameK = '-'.join(str(x) for x in ProdDegRatio) + "K.pdf"
#
# np.savetxt(filenameL, outputMatrixL, delimiter=",")
# np.savetxt(filenameK, outputMatrixK, delimiter=",")
absc = np.arange(0,tmax+dt,dt)

plt.figure()
for l in outputMatrixL:
    print(l.shape, absc.shape)
    plt.plot(absc, l)
plt.legend(ProdDegRatio, loc='upper left')
plt.suptitle("Length of the stem")
# plt.savefig(filenameL, dpi=900)

plt.figure()
for k in outputMatrixK:
    print(l.shape, absc.shape)
    plt.plot(absc, k)
plt.legend(ProdDegRatio, loc='upper left')
plt.suptitle("Curvature of the stem")
# plt.savefig(filenameK, dpi=900)

plt.show()
