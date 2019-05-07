# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 08:24:39 2018

@author: User
"""
#%%

import pandas as pd
import numpy as np
import b_value
import obspy

#%%

april1_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt1\var_april.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
#april1_wfmeas = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt1\var_april.wfmeas', header=None, delim_whitespace=True, usecols=[0,1,8,12], names=['sta', 'chan', 'amp', 'arid'])
#april1_arrival = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt1\var_april.arrival', delim_whitespace=True, header=None, usecols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23], names=['time','arid','jdate','stassid','chanid','chan','iphase','stype','deltim','azimuth','delaz','slow','delslo','ema','rect','amp','per','logat','clip','fm','snr','qual','auth'])

april2_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt2\var_april_after_12.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])

mar_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\march\var_march.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])

feb_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\feb\feb99work.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])

jan_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\jan\pf299.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])

#%%

cat = 