#coding: utf-8
# premier test de modele numerique pour l'instant de nutation par mechanoperception 
# discussion avec Stéphane; as t on besoin d'hysteresis dan sla loi lpoint de sigma ? lag necessaire ? pas forcement? on peut ocmmencer avec juste un truc affine regarder les points de fonctionnement etc ...

import os
import matplotlib.pyplot as plt

import numpy as np


from scipy.optimize import fsolve
import math

#~ def equations(p):
    #~ x, y = p
    #~ return (x+y**2-4, math.exp(x) + x*y - 3)

#~ x, y =  fsolve(equations, (1, 1))

#~ print equations((x, y))




################################################################################################################


#Voir cahier Clairfontaine bleu ciel "Carambola" pour les deux equations d'equilibre mecanique : 
#l = l'ancien ltilde c est la longueur  au milieu d equilibre qu 'on cherche
#k ) kappa la courbure

#~ l10= sont les target de la loi de croissance
#~ l20=






l10=12.
l20=10.
a=0.1
YM=1.0 # pour linstant Young Modulus


alpha=0.001 #a definir 
beta=0.1
tau=5000.0

g0=0.0000 # nearly no basic growth
N1=1.0
N2=1.0

listek=[]
listesigma1=[]
listesigma2=[]
listel10=[]
listel20=[]
listeN1=[]
listeN2=[]
listet=[]
listel=[]

def equations(p):
    l, k = p
    return (3*a*l*(l*(-a*k + 1) - l20)**2 + 3*a*l*(l*(a*k + 1) - l10)**2, -3*(l - l10)**2 + 3*(l - l20)**2 - (-3*a*k + 3)*(l*(-a*k + 1) - l20)**2 + (3*a*k + 3)*(l*(a*k + 1) - l10)**2   )

#~ l, k =  fsolve(equations, (l0, 1.0/20))
#~ print l,k

#~ print equations((l, k))


def stress12(l,l10,l20,a,k):
	#~ sigma1=YM*(l/l10*(1+a*k/2) -1   )
	#~ sigma2=YM*(l/l20*(1-a*k/2) -1   )
	return sigma1,sigma2
	
t=0.0
dt=0.01
tmax=5000
lprevious=0.5*(l10+l20) # initial guess, variable sert always pour initialguess
kprevious=1.0/20 # intial guess
while t <tmax:
	t=t+dt
	# 1. calcul de la nouvelle forme
	l, k =  fsolve(equations, (lprevious, kprevious))
	lprevious=l
	kprevious=k
	
	# 2. calcul des nouveaux stress
	sigma1,sigma2=stress12(l,l10,l20,a,k)

	# update growth rate based on stress : (depend de l'hypothese à tester!!)

	#~ N1=N1*(1 - dt/tau) +max(dt*beta*sigma1,0)
	#~ N2=N2*(1 - dt/tau) +max(dt*beta*sigma2,0)
	
	#~ g1=g0+alpha*N1
	#~ g2=g0+alpha*N2
	
	# grow : 
	l10= l10*(1+max(dt*g1,0))
	l20= l20*(1+max(dt*g2,0))

	listek.append(k)
	listel.append(l)

	listesigma1.append(sigma1)
	listesigma2.append(sigma2)
	listel10.append(l10)
	listel20.append(l20)
	listeN1.append(N1)
	listeN2.append(N2)
	listet.append(t)


plt.figure()
plt.plot(listet,listek)
plt.xlabel('time')
plt.ylabel('curvature')

import pickle
pickle.dump([listet,listek],open("listetlistek.pickle","w"))

plt.figure()
plt.plot(listet,listel10,label='l10')
plt.plot(listet,listel20,label='l20')
plt.plot(np.array(listet),np.array(listel)*(1+np.array(listek)*a),label='l1')
plt.plot(np.array(listet),np.array(listel)*(1-np.array(listek)*a),label='l2')

plt.xlabel('time')
plt.ylabel('lengths')
plt.legend(loc=0)

plt.figure()
plt.plot(listet,listesigma1,label='sigma1')
plt.plot(listet,listesigma2,label='sigma2')
plt.xlabel('time')
plt.ylabel('stress')
plt.legend(loc=0)

plt.figure()
plt.plot(listet,listeN1,label='N1')
plt.plot(listet,listeN2,label='N2')
plt.xlabel('time')
plt.ylabel('N')
plt.legend(loc=0)


plt.show()



