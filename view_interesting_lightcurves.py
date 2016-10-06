import data
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import kicanalysis
"""
kics = [  9592855,   3735597,   8257903,   9101400,   8823397,
   9953894,   9479460,   9357275,   9906590,
   7097571,   8262223,   8043961,   9851944,   6863161,
   6863229]
ka = kicanalysis.kic_analyze(kics)
plt.figure(figsize=(20,16))
for i, kic in enumerate(kics):
    plt.subplot(3,5,i+1)
    ka.lightcurve(kic)
plt.show()
"""
"""
9596089: ~1/3 orbital period
8164262: much faster than ~175 day period
3953106: almost definitely spots (little data)
9479460: ellipsoidal variation? (1/2 orbital period)
6387887: large cool star + hot small?  (same as orbital period)
10028352: likely spots (very small eclipses)
11347875: cepheid of same period as orbit????
"""
kics= [10028352, 6387887, 11347875]
ka = kicanalysis.kic_analyze(kics)
#plt.figure(figsize=(20,16))
for i, kic in enumerate(kics):
    #plt.subplot(3,2,i+1)
    plt.figure(figsize=(14,6))
    plt.subplot(1,2,1)
    ka.lightcurve(kic)
    plt.subplot(1,2,2)
    ka.plot_periodogram(kic)
    plt.show()
