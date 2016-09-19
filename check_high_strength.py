import data
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import kicanalysis

allkics = data.select_kics()

ka = kicanalysis.kic_analyze(allkics)
potential_periods = ka.read_multi_period('per_multi_periods_all_interp_0.7')
array = np.empty((len(potential_periods)*5, 4))
index = 0
for kic, values in potential_periods.iteritems():
    array[index] = [kic, values[0][0], values[0][1], values[5]]
    array[index+1] = [kic, values[1][0], values[1][1], values[5]]
    array[index+2] = [kic, values[2][0], values[2][1], values[5]]
    array[index+3] = [kic, values[3][0], values[3][1], values[5]]
    array[index+4] = [kic, values[4][0], values[4][1], values[5]]
    index +=5
value_array = array.T
kic = value_array[0]
strengths = value_array[2]
mask = np.where(strengths > 0.9)
for No in kic[mask]:
    ka.plot_periodogram(No)