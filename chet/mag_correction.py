# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 13:03:51 2018

@author: cgoerzen
"""
#%%

import obspy
import numpy as np


#%%

def mag_correction(input_ml, sta_lat, sta_lon, ev_lat, ev_lon, depth, s):
    '''
    Implements the magnitude correction as described by Mahani et al., 2018
    
    Params
    ------
    input_ml: The uncorrected magnitude
    sta_lat: The latitude of the receiving station
    sta_lon: The longitude of the receiving station
    ev_lat: The latitude of the earthquake focus
    ev_lon: The longitude of the earthquake focus
    depth: The depth of the earthquake focus
    s: The station correction factor
    
    Returns
    -------
    output_ml: The corrected station magnitude
    '''
    epi_dist = obspy.geodetics.base.gps2dist_azimuth(sta_lat, sta_lon, ev_lat, ev_lon)
    d_hypo = np.sqrt((epi_dist**2) + depth**2)
    if (d_hypo<=85):
        output_ml = [(0.7974*np.log10(d_hypo/100)) + (0.0016*(d_hypo - 100)) + (+2.28566774*np.log10(epi_dist)-1.53265) + 3 + s + input_ml]
    else:
        output_ml = [(-0.1385*np.log10(d_hypo/100)) + (0.0016*(d_hypo - 100)) + (+2.28566774*np.log10(epi_dist)-1.53265) + 3 + s + input_ml]
        
    return output_ml
    
#%%