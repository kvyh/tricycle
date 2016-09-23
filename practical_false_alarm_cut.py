import matplotlib.pyplot as plt
import numpy as np
import periodicity2
import data
import kicanalysis

allkics = data.select_kics()
ka = kicanalysis.kic_analyze(allkics[:1])
ka.interpolate = False
time, x, flux_err, p_orb = ka.get_info(ka.kics[0])
times = time
fluxerr = flux_err

strengths = []
sigma = .001
x = 0
while x < 1000:
    model_data = np.random.normal(1, sigma, len(times))
    pr = periodicity2.Periodicity(times, model_data, fluxerr)
    pr.Lomb_scargle()
    x+=1
    strengths.append(pr.periodogram_results[1])

max_strength = max(strengths)
f = open('practical_false_alarm_cut_constant_mean.txt', 'w')
f.write('Using periodograms with data points cut out periodically (at period '+str(p_orb)+') \n')
f.write('Over the 1000 random-gaussian "lightcurves" used, the strongest period returned by our method in any of them had a strength of '+str(max_strength))

#xs = np.linspace(0, 1, 100001)
#ys = np.array([sum(strengths > x) for x in xs])
#plt.plot(xs, ys)
#plt.xlim((0,.001))
#plt.show()
