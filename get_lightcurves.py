# -*- coding: utf-8 -*-
"""
Select eclipsing binaries from Villanova catalog based on system parameters
and get their light curves.
"""
import glob
import pandas as pd

def select_kics(catFile='villanova-db.csv', dataDir='./data', 
                pmin=0.0, pmax=None):
    '''Return KIC IDs based on system parameters.
    
    Inputs
    ----------
    catFile : string (default: villanova-db.csv)
        name of catalog file
    dataDir : string (default: ./data)
        name of list of light curve filenames
        '''

    '''Construct the query string'''
    qstring = 'period > %s & ' % str(pmin)
    if pmax: qstring += 'period < %s & ' % str(pmax)

    '''Load the Villanova catalog, find the matching KIC IDs.'''
    df = pd.read_csv(catFile)
    kics = df.query(qstring[:-2])['KIC'].values.astype(int).astype(str)

    '''Find light curves for targets.'''   
    filenames = glob.glob(dataDir+'/*.fits')
    lc_list = []
    kics_with_data = []
    for kic in kics:
        files = [f for f in filenames if kic in f]
        lc_list.append(files)
        if len(files) > 0: kics_with_data.append(kic)

    Nfiles = len([item for sublist in lc_list for item in sublist])
    Nsys = len(kics_with_data)

    print 'Found %d systems in catalog meeting criteria, %d of which have data.' % (kics.size, Nsys)    
    print 'Found %d light curves.' % Nfiles

'''Run a test'''    
select_kics(pmin=2.0, pmax=20.)