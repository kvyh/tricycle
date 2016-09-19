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

kics= [7938468, 11347875, 4245897, 3953106, 
5266937, 10028352]
ka = kicanalysis.kic_analyze(kics)
plt.figure(figsize=(20,16))
for i, kic in enumerate(kics):
    plt.subplot(3,2,i+1)
    ka.lightcurve(kic)
plt.show()
