#SINGLE BUNDLE
#Nadrowski - Nick Senofsky

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
Kes = 140 * 10**-6 #stiffness of extent spring
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
#fmax = 439 * 10**-12  #330-800 pN  (352)
#P0 =.2
#P1 = -.8   #"assume that increased Ca levels at the motor site reduce active force generation"
Cmax = 250 * 10**-3  #assuming max calcium concentration @ motor is endolymph concentration
Cm = Cmax * 1

a0 = 1.5
fmax = ((87 * a0) + 352) * 10**-12
P0 = .2
#P1 = -1   #"assume that increased Ca levels at the motor site reduce active force generation"
sp1M = (((0.3)/(87* 10**-12))*(fmax - (352* 10**-12))) + 0.65
P1 = (-1 * P0 * sp1M)/Cm
S = -1 * Cm * P1 / P0
force = fmax / (N_a * P0)

print(S)

#time
tstart = 0
tend = 10
time_N = 500000
dt = (tend - tstart)/time_N
t = np.linspace(tstart,tend,time_N)

#noise
n = 0#math.sqrt(KbT * lambda_ / dt)
na = 0#math.sqrt(Kb * T_a * lambda_a / dt)
dc = 0#.5 * Cm * math.sqrt(tau_c / (N * dt))

z0 = [-5*10**-8, -1.3*10**-7 , .15]



def model(z,time):
    x = z[0]
    xa = z[1]
    c = z[2]

    p0 = (1 + (A * math.exp(-1 * (x - xa) / delta)))**-1
    p = P0 + (P1 * c)

    #External force
    fext = 0 #.52 *10**-9

    Y = x - xa - (D * p0)

    dxdt = a0 * ((1 / lambda_) * ((-1 * Kgs * Y)
                                        - (Ksp * x) + fext + Kgs*(1 + (A * math.exp(0)))**-1))
    dxadt = a0 * ((1 / lambda_a) * ((Kgs * Y)
                                          - (gamma * N_a * force * p) + (Kes * (xa + xes)) + Kgs*(1 + (A * math.exp(0)))**-1
                                    + (gamma * N_a * force * p0)))
    dcdt = a0 * ((1 / tau) * (C0 - c + (Cm * p0)))
    return [dxdt, dxadt, dcdt]

#RK4
def model_RK4():
    X = np.zeros(time_N)
    Xa = np.zeros(time_N)
    C = np.zeros(time_N)
    X[0] = z0[0]
    Xa[0] = z0[1]
    C[0] = z0[2]
    for i in range(0,(time_N - 1)):
        zk1 = model([X[i], Xa[i], C[i]], t[i])
        zk2 = model([X[i] + zk1[0]*dt/2,
                     Xa[i] + zk1[1]*dt/2,
                     C[i] + zk1[2]*dt/2], t[i] + dt/2)
        zk3 = model([X[i] + zk2[0]*dt/2,
                     Xa[i] + zk2[1]*dt/2,
                     C[i] + zk2[2]*dt/2], t[i] + dt/2)
        zk4 = model([X[i] + zk3[0]*dt,
                     Xa[i] + zk3[1]*dt,
                     C[i] + zk3[2]*dt], t[i] + dt)
        X[i + 1] = X[i] + (dt/6)*(zk1[0] + 2*zk2[0]
                                + 2*zk3[0] + zk4[0]) + (dt*n/lambda_)*random.gauss(0, 1.414)
        Xa[i + 1] = Xa[i] + (dt/6)*(zk1[1] + 2*zk2[1]
                                  + 2*zk3[1] + zk4[1]) + (dt*na/lambda_a)*random.gauss(0, 1.414)
        C[i + 1] = C[i] + (dt/6)*(zk1[2] + 2*zk2[2]
                                + 2*zk3[2] + zk4[2]) + (dt*dc/tau)*random.gauss(0, 1.414)
    return [X, Xa, C]


#z = odeint(model, z0, t)
#X = z[:,0]
#Xa = z[:,1]
#C = z[:,2]

z = model_RK4()
X = z[0]
Y = z[1]
TT = z[2]

plt.plot(t,X, "r")
plt.plot(t, Y, "b")

plt.xlabel('Time (s)')
plt.ylabel('Bundle & Motor Displacement (m)')
plt.show()