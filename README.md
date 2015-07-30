3 Periods, 2 Stars, 1 Age

Looking through the *Kepler* [eclipsing binaries](http://keplerebs.villanova.edu), can we find systems with starspot signals from both stars? To find these, look for 3 distinct, non-harmonic, periods.

This would give rotation periods for both components (challenge goal: figure out which period comes from which star).

Follow-up spectroscopic observation then gives masses and radii.

## Retrieving Light Curves
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
