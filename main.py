import os
import numpy as np
from scipy.optimize import fsolve
import math

#================
# Simulation parameters
tmax = 100
dt = 0.01

#================
# Initial parameters
l0 = [10,11,12,13]


#================
# Output Variables
l = np.zeros(shape=(int(tmax/dt)))
