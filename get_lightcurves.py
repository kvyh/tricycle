# -*- coding: utf-8 -*-
"""
Select eclipsing binaries from Villanova catalog based on system parameters
and get their light curves.
"""

def select_kics(catFile='villanova-db.csv', dataDir='./data', 
                pmin=0.0, pmax=None):
    '''
    Return KIC IDs based on system parameters.
    
    Inputs
    ----------
    catFile : string (default: villanova-db.csv)
        name of catalog file
    pmin : float (default: 0.0)
        minimum orbital period in days
    pmax : float (default: None)
        maximum orbital period in days
    '''
    import pandas as pd
    
    '''Construct the query string'''
    qstring = 'period > %s & ' % str(pmin)
    if pmax: qstring += 'period < %s & ' % str(pmax)

    '''Load the Villanova catalog, find the matching KIC IDs.'''
    df = pd.read_csv(catFile)
    kics = df.query(qstring[:-2])['KIC'].values.astype(int).astype(str)

    print 'Found %d systems in catalog meeting criteria.' % kics.size    
    return kics.astype(int)
