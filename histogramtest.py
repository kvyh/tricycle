import data
import matplotlib.pyplot as plt
import numpy as np
import kicanalysis

allkics = data.select_kics()

ka = kicanalysis.kic_analyze(allkics[:])
names = ['periodogram_periods_all_cut_0.4', 'periodogram_periods_all_interp_0.4', 'autocor_periods_all_cut_0.4', 'autocor_periods_all_interp_0.4']
dict1 = ka.readfile(names[0])
dict2 = ka.readfile(names[1])
dict3 = ka.readfile(names[2])
dict4 = ka.readfile(names[3])
print len(dict1)
#ka.histogram(np.array([[dict1,dict2,dict3,dict4],names]))
ka.histogram2D(np.array([[dict1,dict2,dict3,dict4],names]))

