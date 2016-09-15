import data
import matplotlib.pyplot as plt
import numpy as np
import kicanalysis

allkics = data.select_kics()

ka = kicanalysis.kic_analyze(allkics[:])
folder = 'periodicity_results/'
base_names = ['periodogram_periods_all_interp_', 'autocor_periods_all_interp_']
ax = plt.figure(1, (16, 20))
i = 1
for cut in ['0.3','0.4','0.5','0.6','0.65','0.7']:
    files = [folder+base_names[0]+cut,folder+base_names[1]+cut]
    #dict1 = ka.readfile(folder+base_names[0]+cut)
    #dict2 = ka.readfile(folder+base_names[1]+cut)
    better = ka.better_periods(files)
    print len(better)
    plt.subplot(int(str(32)+str(i)))
    ka.plot_results(dictionary=better)
    ka.write(better, 'better', '_interp_'+cut)
    i += 1
plt.show()