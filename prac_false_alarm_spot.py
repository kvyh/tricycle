import matplotlib.pyplot as plt
import numpy as np
import periodicity2
import data
import kicanalysis

allkics = data.select_kics()
ka = kicanalysis.kic_analyze(allkics[:1])
ka.interpolate = True
time, x, flux_err, p_orb = ka.get_info(ka.kics[0])
times = time
fluxerr = flux_err

strengths = {}
sigma = .01
period_orb = [.1, 1,5,10,40,80]
period_rot = [.1,1,2,5,6,10,11,40,42,80,85]

for sigma in sigmas:
    x = 0
    while x < 1000:
        model_data = np.random.normal(1, sigma, len(times))
        pr = periodicity2.Periodicity(times, model_data, fluxerr)
        pr.Lomb_scargle()
        x+=1
        if sigma in strengths:
            strengths[sigma].append(pr.periodogram_results[1])
        else:
            strengths[sigma] = [pr.periodogram_results[1]]

xs = np.linspace(0, 1, 100001)

    ys = np.array([sum(strengths[sigma] > x) for x in xs])
    plt.plot(xs, ys)
plt.xlim((0,.01))
plt.show()

max_strength = max([max(strengths[.001]),max(strengths[1.])])
f = open('practical_false_alarm_spot_model','w')
f.write(')
