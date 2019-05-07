# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 19:50:13 2018

@author: User
"""

#%%

import pandas as pd
import numpy as np
#import vincenty

#%%

richter_cf = pd.read_csv(r'C:\Users\User\Desktop\Work Stuff\Scripts\data\richter_cf.csv')

d_km = richter_cf.d_km
cf = richter_cf.cf

#%%

np.polyfit(np.log10(d_km), cf, 1)

#%%

print(-2.28566774*np.log10(100) + 1.5326484)