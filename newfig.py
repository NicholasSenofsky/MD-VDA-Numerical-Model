#newfig


import numpy as np
import math
from scipy import *
import random
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.gridspec as gs
from scipy.integrate import odeint
from scipy.fft import fft, fftfreq
#time
tstart = 0
tend = 10
time_N = 500000
#time_N = 10000000
dt = (tend - tstart)/time_N
timespace = np.linspace(tstart,tend,time_N)
samplerate = 1/dt
duration = tend - tstart
nb = 15

k15 = np.load(r"D:\HairBundles\4221K15GaussM1P5SP25XLSuprMatrix.npy")
k10 = np.load(r"D:\HairBundles\4221K10GaussM1P5SP25XLSuprMatrix.npy")
k5 = np.load(r"D:\HairBundles\4221K5GaussM1P5SP25XLSuprMatrix.npy")
k = np.load(r"D:\HairBundles\4221NoDamageGaussM1P5SP25XLSuprMatrix.npy")



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

halfsecspace = np.linspace(0,.5,25000)
lessspace = np.linspace(0,.3,15000)

#plt.plot(   linspace(0, 3, len(x)),  x   )  #How to set t0 = 0

#Create Matricies of Open Prob
#p0 = (1 + (A * math.exp(-1 * (x - xa) / delta)))**-1

p=np.zeros(time_N)
for i in range(time_N):
    p[i] = (1 + (A * math.exp(-1 * (k[i,7,7,0] - k[i,7,7,1]) / delta)))**-1

p5=np.zeros(time_N)
for i in range(time_N):
    p5[i] = (1 + (A * math.exp(-1 * (k5[i,7,7,0] - k5[i,7,7,1]) / delta)))**-1

p10=np.zeros(time_N)
for i in range(time_N):
    p10[i] = (1 + (A * math.exp(-1 * (k10[i,7,7,0] - k10[i,7,7,1]) / delta)))**-1

p15=np.zeros(time_N)
for i in range(time_N):
    p15[i] = (1 + (A * math.exp(-1 * (k15[i,7,7,0] - k15[i,7,7,1]) / delta)))**-1

#plots
plt.figure(figsize=(8,8))
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, hspace=0.3, wspace=0.3)
ax1 = plt.subplot2grid((2,2), (0,0), colspan=1, rowspan=1)
ax2 = plt.subplot2grid((2,2), (0,1), colspan=1, rowspan=1)
ax3 = plt.subplot2grid((2,2), (1,0), colspan=1, rowspan=1)
ax4 = plt.subplot2grid((2,2), (1,1), colspan=1, rowspan=1)

ax1.set_xlim([0,.3])
ax2.set_xlim([0,.3])
ax3.set_xlim([0,.3])
ax4.set_xlim([0,.3])
ax1.set_ylim([0,1])
ax2.set_ylim([0,1])
ax3.set_ylim([0,1])
ax4.set_ylim([0,1])

ax1.plot(lessspace, p[185000:200000], "k")
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Transduction Channel Open Probability')

ax2.plot(lessspace, p5[185000:200000], "k")
ax2.set_xlabel('Time (s)')
#ax2.set_ylabel('Transduction Channel Open Probability')

ax3.plot(lessspace, p10[185000:200000], "k")
ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Transduction Channel Open Probability')

ax4.plot(lessspace, p15[185000:200000], "k")
ax4.set_xlabel('Time (s)')
#ax4.set_ylabel('Transduction Channel Open Probability')


ax1.set_yticks([])
ax1.set_xticks([])
ax2.set_yticks([])
ax2.set_xticks([])
ax3.set_yticks([])
ax3.set_xticks([])
ax4.set_yticks([])
ax4.set_xticks([])

scalebar1x = AnchoredSizeBar(ax1.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar2x = AnchoredSizeBar(ax2.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar3x = AnchoredSizeBar(ax3.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar4x = AnchoredSizeBar(ax4.transData, .1, '.1 s', 'lower center', frameon=False,)

ax1.add_artist(scalebar1x)
ax2.add_artist(scalebar2x)
ax3.add_artist(scalebar3x)
ax4.add_artist(scalebar4x)

scaletexty1 = 35500

# scalebar1y = AnchoredSizeBar(ax1.transData, 0.03, '', 'center left', size_vertical=10000, borderpad=2, frameon=False)
# ax1.text(-69, scaletexty1, '10000 timesteps', fontsize=10, va='top',backgroundcolor='w', rotation=90)
# ax1.add_artist(scalebar1y)
#
# scalebar2y = AnchoredSizeBar(ax2.transData, 0.03, '', 'center left', size_vertical=10000, borderpad=2, frameon=False)
# #ax2.text(-69, scaletexty1, '10000 timesteps', fontsize=10, va='top',backgroundcolor='w', rotation=90)
# #ax2.add_artist(scalebar2y)
#
# scalebar3y = AnchoredSizeBar(ax3.transData, 0.03, '', 'center left', size_vertical=10000, borderpad=2, frameon=False)
# ax3.text(-69, scaletexty1, '10000 timesteps', fontsize=10, va='top',backgroundcolor='w', rotation=90)
# ax3.add_artist(scalebar3y)
#
# scalebar4y = AnchoredSizeBar(ax4.transData, 0.03, '', 'center left', size_vertical=10000, borderpad=2, frameon=False)
# #ax4.text(-69, scaletexty1, '10000 timesteps', fontsize=10, va='top',backgroundcolor='w', rotation=90)
# #ax4.add_artist(scalebar4y)

rs1 = ax1.spines["right"]
rs1.set_visible(False)
t1  = ax1.spines["top"]
t1.set_visible(False)
rs2 = ax2.spines["right"]
rs2.set_visible(False)
t2  = ax2.spines["top"]
t2.set_visible(False)
rs3 = ax3.spines["right"]
rs3.set_visible(False)
t3  = ax3.spines["top"]
t3.set_visible(False)
rs4 = ax4.spines["right"]
rs4.set_visible(False)
t4  = ax4.spines["top"]
t4.set_visible(False)


ax1.text(-67, 47000, 'a', fontsize=16, fontweight='bold', va='top')
ax2.text(-67, 47000, 'b', fontsize=16, fontweight='bold', va='top')
ax3.text(-67, 47000, 'c', fontsize=16, fontweight='bold', va='top')
ax4.text(-67, 47000, 'd', fontsize=16, fontweight='bold', va='top')

print("Pmax (k=1) = ")
print(max(p[185000:200000]))
print("Pmax (k=5) = ")
print(max(p5[185000:200000]))
print("Pmax (k=10) = ")
print(max(p10[185000:200000]))
print("Pmax (k=15) = ")
print(max(p15[185000:200000]))


plt.show()
plt.savefig(r'C:\Users\Nick Senofsky\Desktop\fig4.eps', dpi=1200)
