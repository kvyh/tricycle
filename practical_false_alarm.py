import matplotlib.pyplot as plt
import numpy as np
import periodicity2
import data
import kicanalysis

allkics = data.select_kics()
ka = kicanalysis.kic_analyze(allkics[:5])
ka.interpolate = True
for kic in ka.kics:
    time, x, flux_err, p_orb = ka.get_info(kic)
    if len(time) > 60000:
        times = time[:60000]
        fluxerr = flux_err[:60000]

strengths = {}
sigmas = [.001, 1.]
for sigma in sigmas:
    x = 0
    while x < 1000:
        model_data = np.random.normal(1, sigma, 60000)
        pr = periodicity2.Periodicity(times, model_data, fluxerr)
        pr.Lomb_scargle()
        x+=1
        if sigma in strengths:
            strengths[sigma].append(pr.periodogram_results[1])
        else:
            strengths[sigma] = [pr.periodogram_results[1]]

#xs = np.linspace(0, 1, 100001)
#for sigma in sigmas:
#    ys = np.array([sum(strengths[sigma] > x) for x in xs])
#    plt.plot(xs, ys)
#plt.xlim((0,.01))
#plt.show()

max_strength = max([max(strengths[.001]),max(strengths[1.])])
f = open('practical_false_alarm_constant_mean', 'w')
f.write('over the 2000 random-gaussian "lightcurves" used, the strongest period returned by our method in any of them had a strength of '+str(max_strength))
