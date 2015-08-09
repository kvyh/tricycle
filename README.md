![tricycle logo]
(https://github.com/StellarArmy/tricycle/blob/master/logo.png)

3 Periods, 2 Stars, 1 Age

Looking through the *Kepler* [eclipsing binaries](http://keplerebs.villanova.edu), can we find systems with starspot signals from both stars? To find these, look for 3 distinct, non-harmonic, periods.

This would give rotation periods for both components (challenge goal: figure out which period comes from which star).

Follow-up spectroscopic observation then gives masses and radii.

## Retrieving Light Curves ##
*This requires access to the Kepler MySQL database on the UW network.*

```python
import loadlc_db as ldb
import get_lightcurves as glc
import matplotlib.pyplot as plt

'''Select systems with 10.0 < p_orbit < 10.5'''
kics = glc.select_kics(pmin=10.0, pmax=10.05)

'''This returns two KIC IDs'''

'''Examine the light curves'''
for kic in kics:
    time, flux, fluxerr, cadence, quarter, quality = ldb.loadlc_db(kic)
    plt.plot(time, flux)
    plt.show() 
```

## Simulating Light Curves ##
Start by initializing an EB object, inputing the system parameters.
```python
eb = tw.EB(p_orb=7.07, p_rot1=1.51, p_rot2=0.80, depth=0.5, amp_1=0.03, amp_2=0.02, e_dur=0.4)
``` 

Now make a light curve and plot.
```python
time, flux = eb.make_lc(length=100, sig=0.02)
plt.plot(time, flux)
plt.xlim(0,11)
plt.ylim(0,1.2)
plt.xlabel('Time (days)')
plt.ylabel('Relative Flux')
plt.show()
```

![Example Light Curve]
(https://github.com/StellarArmy/tricycle/blob/master/lightcurve.png)

### Parameters ###
----------
* `p_orb` : float - Orbital period in days    
* `p_rot1` : float - Rotation period of star 1
* `p_rot2` : float - Rotation period of star 2
* `depth` : float - Eclipse depth
* `amp_1` : float - Spot amplitude of star 1
* `amp 2` : Spot amplitude of star 2
* `e_dur` : Eclipse duration in days
* `sc` : Default is for long cadence, set to `True` for short cadence
