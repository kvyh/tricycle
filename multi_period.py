import data
import matplotlib.pyplot as plt
import numpy as np
import kicanalysis

allkics = data.select_kics()

ka = kicanalysis.kic_analyze(allkics[:])
ka.l_s_multi_period()
ka.write(ka.periodogram_multi_results, 'per_multi_', '_interp_0.7')
