#Workspace


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


#ave = np.sum(bundles[20000:,:,:,0])/((time_N-20000)*15*15)
#aveaxis = np.repeat(ave, time_N)

#plt.plot(timespace, bundles[:,1,1,0], "b")
#plt.plot(timespace, bundles[:,1,13,0], "g")
#plt.plot(timespace, bundles[:,13,1,0], "y")
#plt.plot(timespace, bundles[:,13,13,0], "r")
#plt.plot(timespace, bundles[:,7,7,0], "k")
#plt.xlabel('Time (s)')
#plt.ylabel('Bundle Position (m)')
#plt.xlim([7.5,8.2])
#plt.ylim([-7*10**-8,-3*10**-8])
#plt.show()


# def move_avg(x, window):
#     if window%2 != 1:
#         print("window must be odd")
#         window += 1
#     filt_x = np.zeros(len(x), dtype=float)
#     for i in range(int(window/2),  len(x) - int(window/2)):  # filter middle of data
#         filt_x[i] = np.sum(x[i - int(window/2):i + int(window/2) + 1])/window
#     for i in range(int(window/2)):  # filter ends of data
#         filt_x[i] = np.sum(x[0:(i + int(window/2))]) / len(x[0:(i + int(window/2))])
#         filt_x[len(x) - 1 - i] = np.sum(x[(len(x) - i - int(window/2)):len(x)]) / len(x[(len(x) - i - int(window/2)):len(x)])
#     return filt_x
#
# def plot_smooth_PSD(x, framerate, filt_window, num_filt):
#     N = len(x)
#     xf = fft(x)
#     f = np.linspace(0.0, framerate/2, int(N/2))
#     #xff = (2.0/(N*framerate)) * np.abs(xf[0:int(N/2)])**2  #use this to find power spectral density
#     xff = (2.0/N)*np.abs(xf[0:int(N/2)])    #use this to find Fourier transform
#     for i in range(num_filt):
#         xff = move_avg(xff, filt_window)
#     return [f, xff]
#     #plt.plot(f, xff)
#     #plt.xlabel('Frequency (Hz)')
#     #plt.ylabel(r'Spectral Density ($nm^2/Hz$)')
#     #plt.ylabel('Amplitude (nm)')
#
# #freqlist = np.zeros((nb,nb))
# #for i in range(nb):
# #    for j in range(nb):
# #        y = bundles[:,i,j,0] - bundles[:,i,j,0].mean()
# #        x = plot_smooth_PSD(y, samplerate, 11, 3)
# #        arr = x[1]
# #        list = arr.tolist()
# #        maximum = max(list)
# #        index = list.index(maximum)
# #        freqlist[i,j] = x[0][index]
# #        print("Just finished" + str(i) + ", " + str(j))
#
# #np.save(r"C:\Users\Nick Senofsky\Desktop\HairBundles\GaussFreqList.npy", freqlist)
#
#freqlist = np.load(r"GaussFreqList.npy")
#freqlistflat = freqlist.flatten()
# for i in range(len(freqlistflat)):
#    print(freqlistflat[i])
#plt.hist(freqlistflat, bins=20, color='k')
#plt.xlabel('Bundle Frequency (Hz)')
#plt.ylabel('Number of Bundles')
#plt.gca().spines["right"].set_visible(False)
#plt.gca().spines["top"].set_visible(False)
#plt.show()

# plt.hist(bundles[:,7,7,0], bins=100, color="k")
# plt.xlabel('Bundle Position (m)')
# plt.ylabel('Histogram')
# plt.xlim([-7*10**-8,-3*10**-8])
# #
# #Pz = np.zeros(time_N)
# #for t in range(0, (time_N)):
# #    Pz[t] = (1 + (A * math.exp(-1 * (bundles[t,1,1,0] - bundles[t,1,1,1]) / delta)))**-1
#
# #plt.plot(timespace, Pz[:], "b")
# #plt.xlabel('Time (s)')
# #plt.ylabel('Open Probability')
# #plt.xlim([3,3.5])"""
#
# 4HISTOGRAM
# plt.figure(figsize=(8,8))
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, hspace=0.3, wspace=0.3)
# ax1 = plt.subplot2grid((2,2), (0,0), colspan=1, rowspan=1)
# ax2 = plt.subplot2grid((2,2), (0,1), colspan=1, rowspan=1)
# ax3 = plt.subplot2grid((2,2), (1,0), colspan=1, rowspan=1)
# ax4 = plt.subplot2grid((2,2), (1,1), colspan=1, rowspan=1)
#
# ax1.set_xlim([-70,-30])
# ax2.set_xlim([-70,-30])
# ax3.set_xlim([-70,-30])
# ax4.set_xlim([-70,-30])
# ax1.set_ylim([0,50000])
# ax2.set_ylim([0,50000])
# ax3.set_ylim([0,50000])
# ax4.set_ylim([0,50000])
#
# binlist = np.linspace(-70,-30,150)
#
# ax1.hist(k[0:500000,7,7,0]*10**9, binlist, color="k")
# #ax1.set_xlabel('Bundle Position (nm)')
# ax1.set_ylabel('No. of Timesteps')
#
# ax2.hist(k5[0:500000,7,7,0]*10**9,binlist, color="k")
# #ax2.set_xlabel('Bundle Position (nm)')
# #ax2.set_ylabel('No. of Timesteps')
#
# ax3.hist(k10[0:500000,7,7,0]*10**9,binlist, color="k")
# ax3.set_xlabel('Bundle Position (nm)')
# ax3.set_ylabel('No. of Timesteps')
#
# ax4.hist(k15[0:500000,7,7,0]*10**9,binlist, color="k")
# ax4.set_xlabel('Bundle Position (nm)')
# #ax4.set_ylabel('No. of Timesteps')
#
# ax1.set_yticks([])
# ax1.set_xticks([])
# ax2.set_yticks([])
# ax2.set_xticks([])
# ax3.set_yticks([])
# ax3.set_xticks([])
# ax4.set_yticks([])
# ax4.set_xticks([])
#
# scalebar1x = AnchoredSizeBar(ax1.transData, 20, '20 nm', 'upper center', frameon=False,)
# scalebar2x = AnchoredSizeBar(ax2.transData, 20, '20 nm', 'upper center', frameon=False,)
# scalebar3x = AnchoredSizeBar(ax3.transData, 20, '20 nm', 'upper center', frameon=False,)
# scalebar4x = AnchoredSizeBar(ax4.transData, 20, '20 nm', 'upper center', frameon=False,)
#
# ax1.add_artist(scalebar1x)
# ax2.add_artist(scalebar2x)
# #ax3.add_artist(scalebar3x)
# #ax4.add_artist(scalebar4x)
#
# scaletexty1 = 35500
#
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
#
# rs1 = ax1.spines["right"]
# rs1.set_visible(False)
# t1  = ax1.spines["top"]
# t1.set_visible(False)
# rs2 = ax2.spines["right"]
# rs2.set_visible(False)
# t2  = ax2.spines["top"]
# t2.set_visible(False)
# rs3 = ax3.spines["right"]
# rs3.set_visible(False)
# t3  = ax3.spines["top"]
# t3.set_visible(False)
# rs4 = ax4.spines["right"]
# rs4.set_visible(False)
# t4  = ax4.spines["top"]
# t4.set_visible(False)
#
#
# ax1.text(-67, 47000, 'a', fontsize=16, fontweight='bold', va='top')
# ax2.text(-67, 47000, 'b', fontsize=16, fontweight='bold', va='top')
# ax3.text(-67, 47000, 'c', fontsize=16, fontweight='bold', va='top')
# ax4.text(-67, 47000, 'd', fontsize=16, fontweight='bold', va='top')
#
# plt.show()
#plt.savefig(r'C:\Users\Nick Senofsky\Desktop\fig4.eps', dpi=1200)

TRACES
halfsecspace = np.linspace(0,.5,25000)
lessspace = np.linspace(0,.3,15000)

#plt.plot(   linspace(0, 3, len(x)),  x   )  #How to set t0 = 0

 plt.figure(figsize=(8,20))
 plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, hspace=0.3, wspace=0.3)
 ax1 = plt.subplot2grid((4,2), (0,0), colspan=1, rowspan=1)
 ax2 = plt.subplot2grid((4,2), (1,0), colspan=1, rowspan=1)
 ax3 = plt.subplot2grid((4,2), (2,0), colspan=1, rowspan=1)
 ax4 = plt.subplot2grid((4,2), (3,0), colspan=1, rowspan=1)
 ax5 = plt.subplot2grid((4,2), (0,1), colspan=1, rowspan=1)
 ax6 = plt.subplot2grid((4,2), (1,1), colspan=1, rowspan=1)
 ax7 = plt.subplot2grid((4,2), (2,1), colspan=1, rowspan=1)
 ax8 = plt.subplot2grid((4,2), (3,1), colspan=1, rowspan=1)

 ax1.set_xlim([0,.5])
 ax2.set_xlim([0,.5])
 ax3.set_xlim([0,.5])
 ax4.set_xlim([0,.5])
 ax5.set_xlim([0,.3])
 ax6.set_xlim([0,.3])
 ax7.set_xlim([0,.3])
 ax8.set_xlim([0,.3])
 ax1.set_ylim([-70,30])
 ax2.set_ylim([-70,30])
 ax3.set_ylim([-70,30])
 ax4.set_ylim([-70,30])
 ax5.set_ylim([-70,-30])
 ax6.set_ylim([-70,-30])
 ax7.set_ylim([-70,-30])
 ax8.set_ylim([-70,-30])

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
 rs5 = ax5.spines["right"]
 rs5.set_visible(False)
 t5  = ax5.spines["top"]
 t5.set_visible(False)
 rs6 = ax6.spines["right"]
 rs6.set_visible(False)
 t6  = ax6.spines["top"]
 t6.set_visible(False)
 rs7 = ax7.spines["right"]
 rs7.set_visible(False)
 t7  = ax7.spines["top"]
 t7.set_visible(False)
 rs8 = ax8.spines["right"]
 rs8.set_visible(False)
 t8  = ax8.spines["top"]
 t8.set_visible(False)


 ax1.plot(halfsecspace, k[175000:200000,1,1,0]*10**9 + 60, "lime")
 ax1.plot(halfsecspace, k[175000:200000,1,13,0]*10**9 + 45, "green")
 ax1.plot(halfsecspace, k[175000:200000,13,1,0]*10**9 + 30, "lawngreen")
 ax1.plot(halfsecspace, k[175000:200000,13,13,0]*10**9 + 15, "darkgreen")
 ax1.plot(halfsecspace, k[175000:200000,7,7,0]*10**9, "k")
 ax1.set_xlabel('Time (s)')
 ax1.set_ylabel('Bundle Position (nm)')

 ax5.plot(lessspace, k[185000:200000,7,7,0]*10**9, "k")
 ax5.set_xlabel('Time (s)')
 ax5.set_ylabel('Bundle Position (nm)')

ax2.plot(halfsecspace, k5[375000:400000,1,1,0]*10**9 + 60, "lime")
ax2.plot(halfsecspace, k5[375000:400000,1,13,0]*10**9 + 45, "green")
ax2.plot(halfsecspace, k5[375000:400000,13,1,0]*10**9 + 30, "lawngreen")
ax2.plot(halfsecspace, k5[375000:400000,13,13,0]*10**9 + 15, "darkgreen")
ax2.plot(halfsecspace, k5[375000:400000,7,7,0]*10**9, "k")
#ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Bundle Position (nm)')

ax6.plot(lessspace, k5[185000:200000,7,7,0]*10**9, "k")
#ax6.set_xlabel('Time (s)')
#ax6.set_ylabel('Bundle Position (nm)')

ax3.plot(halfsecspace, k10[175000:200000,1,1,0]*10**9 + 60, "lime")
ax3.plot(halfsecspace, k10[175000:200000,1,13,0]*10**9 + 45, "green")
ax3.plot(halfsecspace, k10[175000:200000,13,1,0]*10**9 + 30, "lawngreen")
ax3.plot(halfsecspace, k10[175000:200000,13,13,0]*10**9 + 15, "darkgreen")
ax3.plot(halfsecspace, k10[175000:200000,7,7,0]*10**9, "k")
#ax3.set_xlabel('Time (s)')
ax3.set_ylabel('Bundle Position (nm)')

ax7.plot(lessspace, k10[185000:200000,7,7,0]*10**9, "k")
#ax7.set_xlabel('Time (s)')
#ax7.set_ylabel('Bundle Position (nm)')

ax4.plot(halfsecspace, k15[175000:200000,1,1,0]*10**9 + 60, "lime")
ax4.plot(halfsecspace, k15[175000:200000,1,13,0]*10**9 + 45, "green")
ax4.plot(halfsecspace, k15[175000:200000,13,1,0]*10**9 + 30, "lawngreen")
ax4.plot(halfsecspace, k15[175000:200000,13,13,0]*10**9 + 15, "darkgreen")
ax4.plot(halfsecspace, k15[175000:200000,7,7,0]*10**9, "k")
ax4.set_xlabel('Time (s)')
ax4.set_ylabel('Bundle Position (nm)')

ax8.plot(lessspace, k15[185000:200000,7,7,0]*10**9, "k")
ax8.set_xlabel('Time (s)')
#ax8.set_ylabel('Bundle Position (nm)')

ax1.text(.03, 28, 'a', fontsize=16, fontweight='bold', va='top')
ax2.text(.03, 28, 'b', fontsize=16, fontweight='bold', va='top')
ax3.text(.03, 28, 'c', fontsize=16, fontweight='bold', va='top')
ax4.text(.03, 28, 'd', fontsize=16, fontweight='bold', va='top')
ax5.text(.02, -31, 'e', fontsize=16, fontweight='bold', va='top')
ax6.text(.02, -31, 'f', fontsize=16, fontweight='bold', va='top')
ax7.text(.02, -31, 'g', fontsize=16, fontweight='bold', va='top')
ax8.text(.02, -31, 'h', fontsize=16, fontweight='bold', va='top')

ax1.set_yticks([])
ax1.set_xticks([])
ax2.set_yticks([])
ax2.set_xticks([])
ax3.set_yticks([])
ax3.set_xticks([])
ax4.set_yticks([])
ax4.set_xticks([])
ax5.set_yticks([])
ax5.set_xticks([])
ax6.set_yticks([])
ax6.set_xticks([])
ax7.set_yticks([])
ax7.set_xticks([])
ax8.set_yticks([])
ax8.set_xticks([])

scalebar1x = AnchoredSizeBar(ax1.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar2x = AnchoredSizeBar(ax2.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar3x = AnchoredSizeBar(ax3.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar4x = AnchoredSizeBar(ax4.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar5x = AnchoredSizeBar(ax5.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar6x = AnchoredSizeBar(ax6.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar7x = AnchoredSizeBar(ax7.transData, .1, '.1 s', 'lower center', frameon=False,)
scalebar8x = AnchoredSizeBar(ax8.transData, .1, '.1 s', 'lower center', frameon=False,)

#ax1.add_artist(scalebar1x)
#ax2.add_artist(scalebar2x)
#ax3.add_artist(scalebar3x)
ax4.add_artist(scalebar4x)
#ax5.add_artist(scalebar5x)
#ax6.add_artist(scalebar6x)
#ax7.add_artist(scalebar7x)
ax8.add_artist(scalebar8x)

scaletexty1 = -6
scaletexty2 = -38

xord1 = .039
point1 = [xord1, -30]
point2 = [xord1, -10]
x_values1 = [point1[0],point2[0]]
y_values1 = [point1[1],point2[1]]

xord2 = .023
point3 = [xord2, -60]
point4 = [xord2, -40]
x_values2 = [point3[0],point4[0]]
y_values2 = [point3[1],point4[1]]


#scalebar1y = AnchoredSizeBar(ax1.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax1.text(.012, scaletexty1, '20 nm', fontsize=10, va='top', rotation=90, c = 'black', backgroundcolor = 'w', zorder = 7)
ax1.plot(x_values1, y_values1, c = 'black', zorder = 8)
#ax1.add_artist(scalebar1y)

#scalebar2y = AnchoredSizeBar(ax2.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax2.text(.012, scaletexty1, '20 nm', fontsize=10, va='top',c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax2.plot(x_values1, y_values1, c = 'black', zorder = 8)
#ax2.add_artist(scalebar2y)

#scalebar3y = AnchoredSizeBar(ax3.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax3.text(.012, scaletexty1, '20 nm', fontsize=10, va='top',c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax3.plot(x_values1, y_values1, c = 'black', zorder = 8)
#ax3.add_artist(scalebar3y)

#scalebar4y = AnchoredSizeBar(ax4.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax4.text(.012, scaletexty1, '20 nm', fontsize=10, va='top',c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax4.plot(x_values1, y_values1, c = 'black', zorder = 8)
#ax4.add_artist(scalebar4y)

#scalebar5y = AnchoredSizeBar(ax5.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax5.text(.0077, scaletexty2, '    20 nm     ', fontsize=10, va='top',c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax5.plot(x_values2, y_values2, c = 'black', zorder = 8)
#ax5.add_artist(scalebar5y)

#scalebar6y = AnchoredSizeBar(ax6.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax6.text(.0077, scaletexty2, '    20 nm     ', fontsize=10, va='top',c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax6.plot(x_values2, y_values2, c = 'black', zorder = 8)
#ax6.add_artist(scalebar6y)

#scalebar7y = AnchoredSizeBar(ax7.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax7.text(.0077, scaletexty2, '    20 nm     ', fontsize=10, va='top', c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax7.plot(x_values2, y_values2, c = 'black', zorder = 8)
#ax7.add_artist(scalebar7y)

#scalebar8y = AnchoredSizeBar(ax8.transData, 0.0001, '', 'center left', size_vertical=20, borderpad=2, frameon=False)
ax8.text(.0077, scaletexty2, '    20 nm     ', fontsize=10, va='top', c = 'black', rotation=90, backgroundcolor = 'w', zorder = 7)
ax8.plot(x_values2, y_values2, c = 'black', zorder = 8)
#ax8.add_artist(scalebar8y)



plt.show()

