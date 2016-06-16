import numpy as np
def curve_cut(self, widthfactor = 0.5):
        """This a function that takes parameters of a star and gives you the time and flux with the eclipses cut out
    
        Parameters:
        Time
        Flux
        bjdo
        period
        pwidth
        swidth
        separation
        
        return:
            timecut and fluxcut"""
        # Convert BJD to Kepler date
        t_0 = self.t_0 - 54833

        phase=((self.time- t_0) % self.p_orb) / self.p_orb
        mask= ((phase > self.p_width * widthfactor) &
               (phase < 1 - self.p_width * widthfactor)) & \
              ((phase > self.sep + self.s_width * widthfactor) |
               (phase < self.sep - self.s_width * widthfactor))
               
        timecut= self.time[mask]
        fluxcut= self.flux[mask]
        errcut = self.err[mask]
        quartercut = self.quarter[mask]

        return timecut, fluxcut, errcut, quartercut
        
def cut_and_intorpolate(self, linear = True):
    '''
    
    '''
    if linear:
        # Convert BJD to Kepler date
        t_0 = self.t_0 - 54833

        phase=((self.time- t_0) % self.p_orb) / self.p_orb

        mask= ((phase > self.p_width * widthfactor) &
               (phase < 1 - self.p_width * widthfactor)) & \
              ((phase > self.sep + self.s_width * widthfactor) |
               (phase < self.sep - self.s_width * widthfactor))
        flux = np.copy(self.flux)    
        flux[~mask] = np.interp(self.time[~mask], self.time[mask], self.flux[mask])

    return self.time, flux, self.err, self.quarter