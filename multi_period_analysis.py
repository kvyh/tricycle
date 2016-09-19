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
periods = value_array[1]
strengths = value_array[2]
p_orb = value_array[3]
p_orb_rot = p_orb/periods
print max(strengths)
plt.figure(figsize=(24, 10))

orb_rot_by_strength = p_orb_rot
p_orb_by_strength = p_orb

plt.subplot(121)
plt.scatter(p_orb_by_strength, orb_rot_by_strength,
            c=strengths, cmap="Blues", s=5, linewidths=0)
plt.title('higher signal strength = darker blue')
plt.ylabel('orbital period/rotation period')
plt.xlabel('orbital period')
plt.xlim((0, 100))
plt.ylim((0, 6))


min_str = 0.15
strength_mask = np.where(min_str <= strengths)
p_orb_rot = p_orb_rot[strength_mask]
p_orb = p_orb[strength_mask]

#ratio_mask = np.where((p_orb_rot < .99) | (p_orb_rot > 1.01))
#p_orb_rot = p_orb_rot[ratio_mask]
#p_orb = p_orb[ratio_mask]

plt.subplot(122)
plt.plot(p_orb, p_orb_rot, linestyle='none', marker='o', markersize=1)
plt.title('target periods with strength greater than ' + str(min_str))
plt.ylabel('orbital period/rotation period')
plt.xlabel('orbital period')
plt.xlim((0, 100))
plt.ylim((0, 6))
plt.show()
