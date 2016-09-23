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
periods = {}
sigma = .1
period_rot = [.1,1,10,40,80]
amplitude = [.01, .05, 1., 3., 10.]
for period in period_rot:
    for amp in amplitude:
        x = 0
        strengths[str(period)+ '_' +str(amp)]=[]
        periods[str(period)+ '_' +str(amp)]=[]
        while x < 100:
            model_data = np.random.normal(1, sigma, len(times)) + amp*np.sin(2*np.pi*times/period)
            pr = periodicity2.Periodicity(times, model_data, fluxerr)
            pr.Lomb_scargle()
            x+=1
            strengths[str(period)+ '_' +str(amp)].append(pr.periodogram_results[1])
            periods[str(period)+ '_' +str(amp)].append(pr.periodogram_results[0])
f = open('practical_false_alarm_sinusoid','w')
f.write('using a sinusoid in addition to gaussian noise (100 cases per period+amplitude combination)\n')
for key, values in strengths.iteritems():
    f.write('for period_amplitude of '+key+', the strongest period was ' +str(periods[key][np.argmax(values)]) + 'with a strength of ' + str(np.max(values)) + '\n')

#xs = np.linspace(0, 1, 100001)
#ys = np.array([sum(strengths[sigma] > x) for x in xs])
#plt.plot(xs, ys)
#plt.xlim((0,.01))
#plt.show()



