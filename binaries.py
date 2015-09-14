"""
This is a module to simulate and analyze light curves of 
spotted eclipsing binaries. 
"""
import numpy as np

class EB(object):
    '''
    Creat an eclipsing binary system object.

    Attributes
    ----------
    p_orb : float
        Orbital period in days    
    p_rot1 : float
        Rotation period of star 1
    p_rot2 : float
        Rotation period of star 2
    depth : float 
        Eclipse depth
    amp_1 : float
        Spot amplitude of star 1
    amp 2 : float
        Spot amplitude of star 2
    e_dur : float
        Eclipse duration in days
    sc : bool
        Default is for long cadence, set to True for short cadence
    '''
        
    def __init__(self, p_orb=7.07, p_rot1=1.51, p_rot2=0.80, depth=0.5, 
                 amp_1=0.03, amp_2=0.02, e_dur=0.4, sc=False):
        self.p_orb = p_orb
        self.p_rot1 = p_rot1
        self.p_rot2 = p_rot2
        self.depth = depth
        self.amp_1 = amp_1
        self.amp_2 = amp_2
        self.e_dur = e_dur
        self.sc = sc

    def make_lc(self, length=100.0, sig=0.02):
        '''
        Make a light curve.
        
        Arguments
        ---------
        length : float
            Length of the light curve in days
        sig : float
            Standard deviation of flux errors
            
        Returns
        -------
        time : array
            Time array
        flux : array
            Flux array
        '''
        
        if self.sc: exptime = 1.0 / 60.0 / 24.0
        else: exptime = 30.0 / 60.0 / 24.0
        time = np.arange(0, 100, exptime)
        t_0 = min(time) + self.p_orb/3. # midpoint of 1st eclipse
        
        '''Make the starspot signals'''
        f_1 = -self.amp_1 * np.sin(time * (2*np.pi) / self.p_rot1)
        f_2 = -self.amp_2 * np.sin(time * (2*np.pi) / self.p_rot2)

        '''Make the eclipse signal as a step function'''
        in_e = np.where((((time - t_0) % self.p_orb)/self.p_orb < (self.e_dur/self.p_orb)))
        f_e = np.zeros_like(time)
        f_e[in_e] = -self.depth
        
        '''Make random noise.'''
        noise = np.random.random(len(time)) * sig
        
        '''Add signals together'''
        flux = 1 + f_e + f_1 + f_2 + noise 
        
        return time, flux