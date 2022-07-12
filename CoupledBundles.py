#CoupledBundles - Nick Senofsky

import numpy as np
import math
from scipy import *
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint


#Variables
lambda_ = 2.8 * 10**-6 #friction coeff of hair bundle
lambda_a = 10 * 10**-6 #friction coeff of adaptation motors
Kgs = 750 * 10**-6 #combinded gating spring stiffness
Ksp = 600 * 10**-6 #combined stiffness of stereociliary pivots and load
Kes = 0 #140 * 10**-6 #stiffness of extent spring
xes = 20 * 10**-9
d = 8.7 * 10**-9 #gating-spring elongation on channel opening
gamma = 0.14 #geometrical gain of stereociliary shear motion
tau = 0.1 * 10**-3 #time constant of calcium feedback
tau_c =  1 * 10**-3 #dwell time of transduction channels
C0 = 0 #intracellular Ca2+ concentration with closed channels
N = 50 ## of stereocillia
N_a = 3000 ## of motors in the hair bundle
Kb = 1.38 * 10**-23 #boltman
T = 300 #temp
KbT = 4.14 * 10**-21
dG = 4.14 * 10**-20 #Gibbs free energy
T_a = 1.5 * T
D = (d / gamma)
A = math.exp((dG + ((Kgs * D**2) / (2*N))) / KbT)
delta = N * KbT / (Kgs * D)
#fmax = 550 * 10**-12  #330-800 pN  (352)
#P0 =.2
#P1 = -1   #"assume that increased Ca levels at the motor site reduce active force generation"
#fmax = 439 * 10**-12  #330-800 pN  (352)
#P0 =.2
#P1 = -.8   #"assume that increased Ca levels at the motor site reduce active force generation"
Cmax = 250 * 10**-3  #assuming max calcium concentration @ motor is endolymph concentration
Cm = Cmax * 1
#S = -1 * Cm * P1 / P0
distance = 50 * 10**-6 #distance between hair bundles 50 um
#K = .003 #(3 pN/nm)
K = .0014 #combined sti↵ness of the membrane, hair bundles, and fillaments
#print(S)

alt = 10
nb = 15
starti = 4
endi = 10

#time
tstart = 0
tend = 10
time_N = 500000
#time_N = 10000000
dt = (tend - tstart)/time_N

#noise
n = math.sqrt(KbT * lambda_ / dt)
na = math.sqrt(Kb * T_a * lambda_a / dt)
dc = 0.5 * Cm * math.sqrt(tau_c / (N * dt))

z0 = [-5*10**-8, -1.3*10**-7 , .15]


bundles = np.zeros((time_N, nb, nb, 3))   #4-d vector to represent our 2 d vector of 1 d vectors through time

b0 = []
for i in range(nb):
    row = []
    for j in range(nb):
        row.append(z0)
    b0.append(row)

bundles[0] = b0

#a0 = np.random.rand(nb,nb)
#for i in range(nb):
#    for j in range(nb):
#        a0[i,j] = 1 #(.5 * a0[i,j]) + 1

a0 = np.random.normal(1.5,.25,size=(nb,nb))
#for i in range(nb):
#    for j in range(nb):
#        a0[i,j] = a0[i,j]

fmaxM = np.random.rand(nb,nb)
for i in range(nb):
    for j in range(nb):
        fmaxM[i,j] = ((87 * a0[i,j]) + 352) * 10**-12

P0 = .2

P1M = np.zeros((nb,nb))
for i in range(nb):
    for j in range(nb):
        sp1M = (((0.3)/(87* 10**-12))*(fmaxM[i,j] - (352* 10**-12))) + 0.65
        P1M[i,j] = (-1 * P0 * sp1M)/Cm

forceM = np.zeros((nb,nb))
for i in range(nb):
    for j in range(nb):
        forceM[i,j] = fmaxM[i,j] / (N_a * P0)


fspring = np.zeros((nb,nb)) #store the spring external forces on each bundle here

def model(z,fkx,t, i, j):
    x = z[0]
    xa = z[1]
    c = z[2]

    p0 = (1 + (A * math.exp(-1 * (x - xa) / delta)))**-1
    p = P0 + (P1M[i,j] * c)

    fext = 0
    #fext = .52 * 10 ** -9
    #if((starti < i < endi)&(starti < j < endi)):
    #    fext= .52*10**-9/alt


    Y = x - xa - (D * p0)

    dxdt = a0[i,j] * ((1/lambda_)*((-1 * Kgs * Y)
                        - (Ksp * x) + fkx +fext))
    dxadt = a0[i,j] * ((1 / lambda_a) * ((Kgs * Y)
                              - (gamma * N_a * forceM[i,j] * p) + (Kes * (xa + xes))))
    dcdt = a0[i,j] * ((1/tau)*(C0 - c + (Cm * p0)))
    return [dxdt, dxadt, dcdt]

#RK4
def model_RK4(z,fkx,t,i,j):
    X = z[0]
    Xa = z[1]
    C = z[2]
    zk1 = model([X, Xa, C], fkx,t,i,j)
    zk2 = model([X + zk1[0]*dt/2,
                Xa + zk1[1]*dt/2,
                C + zk1[2]*dt/2], fkx,t,i,j)
    zk3 = model([X + zk2[0]*dt/2,
                Xa + zk2[1]*dt/2,
                C + zk2[2]*dt/2], fkx,t,i,j)
    zk4 = model([X + zk3[0]*dt,
                Xa + zk3[1]*dt,
                C + zk3[2]*dt], fkx,t,i,j)
    Xn = X + (dt/6)*(zk1[0] + 2*zk2[0]
                                + 2*zk3[0] + zk4[0]) + (dt*n/lambda_)*random.gauss(0, 1.414)
    Xan = Xa + (dt/6)*(zk1[1] + 2*zk2[1]
                                  + 2*zk3[1] + zk4[1]) + (dt*na/lambda_a)*random.gauss(0, 1.414)
    Cn = C + (dt/6)*(zk1[2] + 2*zk2[2]
                                + 2*zk3[2] + zk4[2]) + (dt*dc/tau)*random.gauss(0, 1.414)
    return [Xn, Xan, Cn]

def calcspring(Kparam, k, l, i, j, p, q, t):
    deltaX = bundles[t, p, q, 0] - bundles[t, i, j, 0]
    L0 = math.sqrt(((k**2) + (l**2)) * distance**2)
    F = (Kparam * (1 - (L0/math.sqrt((deltaX + (k*distance))**2 + (l * distance)**2)))
        * (deltaX + (k * distance)))
    return F

def updatefspringdamage(t):
    kalt = K/alt
    for i in range(0, nb):
        for j in range(0, nb):
            spring = 0
            if((i-1 >= 0) and (j-1 >= 0)):
                if((starti<i<endi)and(starti<j<endi)and(starti<(i-1)<endi)and(starti<(j-1)<endi)):
                    spring = spring + calcspring(kalt,-1,1,i,j,i-1,j-1,t)
                else:
                    spring = spring + calcspring(K,-1,1,i,j,i-1,j-1,t)
            if(i-1 >= 0):
                if((starti<i<endi)and(starti<j<endi)and(starti<(i-1)<endi)):
                    spring = spring + calcspring(kalt,-1,0,i,j,i-1,j,t)
                else:
                    spring = spring + calcspring(K,-1,0,i,j,i-1,j,t)
            if((i-1 >= 0) and (j+1 <= nb-1)):
                if((starti<i<endi)and(starti<j<endi)and(starti<(i-1)<endi)and(starti<(j+1)<endi)):
                    spring = spring + calcspring(kalt,-1,1,i,j,i-1,j+1,t)
                else:
                    spring = spring + calcspring(K,-1,1,i,j,i-1,j+1,t)
            if(j-1 >= 0):
                if ((starti<i<endi)and(starti<j<endi)and(starti<(j-1)<endi)):
                    spring = spring + calcspring(kalt,0,1,i,j,i,j-1,t)
                else:
                    spring = spring + calcspring(K,0,1,i,j,i,j-1,t)
            if(j+1 <= nb-1):
                if ((starti<i<endi)and(starti<j<endi)and(starti<(j+1)<endi)):
                    spring = spring + calcspring(kalt,0,1,i,j,i,j+1,t)
                else:
                    spring = spring + calcspring(K,0,1,i,j,i,j+1,t)
            if((i+1 <= nb-1) and (j-1 >= 0)):
                if ((starti<i<endi)and(starti<j<endi)and(starti<(i+1)<endi)and(starti<(j-1)<endi)):
                    spring = spring + calcspring(kalt,1,1,i,j,i+1,j-1,t)
                else:
                    spring = spring + calcspring(K,1,1,i,j,i+1,j-1,t)
            if(i+1 <= nb-1):
                if ((starti<i<endi)and(starti<j<endi)and(starti<(i+1)<endi)):
                    spring = spring + calcspring(kalt,1,0,i,j,i+1,j,t)
                else:
                    spring = spring + calcspring(K,1,0,i,j,i+1,j,t)
            if((i+1 <= nb-1) and (j+1 <= nb-1)):
                if ((starti<i<endi)and(starti<j< endi)and(starti<(i+1)<endi)and(starti<(j+1)<endi)):
                    spring = spring + calcspring(kalt,1,1,i,j,i+1,j+1,t)
                else:
                    spring = spring + calcspring(K,1,1,i,j,i+1,j+1,t)
            fspring[i,j] = spring



def cross_correlation(x, y):
    assert len(x) == len(y)
    X = x - np.sum(x)/len(x)
    Y = y - np.sum(y)/len(y)
    X /= np.std(X)
    Y /= np.std(Y)
    return np.sum(X*Y)/len(X)


def simulation():
    for t in range(0, (time_N -1)):
        updatefspringdamage(t)
        if (t==(time_N/4)):
            print("One Quarter Done")
        if(t==(time_N/2)):
            print("Half way there")
        if (t==(3*time_N/4)):
            print("Just a lil longer")
        for i in range(0, nb):
            for j in range(0, nb):
                bundles[t+1,i,j] = model_RK4(bundles[t,i,j],fspring[i,j],t,i,j)

timespace = np.linspace(tstart,tend,time_N)

def correlate():
    H1H2 = cross_correlation(bundles[:, 1, 1, 0], bundles[:, 1, 5, 0])
    H1H3 = cross_correlation(bundles[:, 1, 1, 0], bundles[:, 5, 1, 0])
    H1H4 = cross_correlation(bundles[:, 1, 1, 0], bundles[:, 5, 5, 0])
    H2H3 = cross_correlation(bundles[:, 1, 5, 0], bundles[:, 5, 1, 0])
    H2H4 = cross_correlation(bundles[:, 1, 5, 0], bundles[:, 5, 5, 0])
    H3H4 = cross_correlation(bundles[:, 5, 1, 0], bundles[:, 5, 5, 0])
    HAC = (H1H2 + H1H3 + H1H4 + H2H3 + H2H4 + H3H4)/6
    D1 = cross_correlation(bundles[:, 3, 3, 0], bundles[:, 1, 1, 0])
    D2 = cross_correlation(bundles[:, 3, 3, 0], bundles[:, 1, 5, 0])
    D3 = cross_correlation(bundles[:, 3, 3, 0], bundles[:, 5, 1, 0])
    D4 = cross_correlation(bundles[:, 3, 3, 0], bundles[:, 5, 5, 0])
    DAC = (D1 + D2 + D3 + D4)/4
    print('The average healthy cross correlation is:')
    print(HAC)
    print('The average damaged cross correlation is:')
    print(DAC)

def outofbounds():
    #reaches steady state in about .5 s
    avstdh = (np.std(bundles[2500:,1,1,0]) + np.std(bundles[2500:,1,5,0]) +
                np.std(bundles[2500:, 5, 1, 0]) + np.std(bundles[2500:,5,5,0]))/4
    avh = (np.mean(bundles[2500:,1,1,0]) + np.mean(bundles[2500:,1,5,0]) +
                np.mean(bundles[2500:, 5, 1, 0]) + np.mean(bundles[2500:,5,5,0]))/4
    dob = 0
    for t in range(2500, (time_N -1)):
        if ((bundles[t,3,3,0] < (avh - 4*avstdh)) or (bundles[t,3,3,0] > (avh + 4*avstdh))):
            dob = dob + 1
    print("Damaged out of bounds for this many seconds:")
    dobs = dob / 50000
    print(dobs)


print("numerical model calculating")

simulation()

#correlate()
#outofbounds()

print("numerical model complete")

print("Saving")

np.save(r"D:\HairBundles\4221K10GaussM1P5SP25NoKesXLSuprMatrix.npy", bundles)

#plt.plot(timespace, bundles[:,3,3,0], "k")
#plt.plot(timespace, bundles[:,1,1,0], "b")
#plt.plot(timespace, bundles[:,1,5,0], "g")
#plt.plot(timespace, bundles[:,5,1,0], "y")
#plt.plot(timespace, bundles[:,5,5,0], "m")
#plt.xlabel('Time (s)')
#plt.ylabel('Bundle Displacement (m)')
#plt.show()

print("Done Saving")

