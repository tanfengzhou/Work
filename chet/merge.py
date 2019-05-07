#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:46:56 2019

@author: pgcseismolab
"""

#%%

import pandas as pd

#%%

dat_path = '~/Desktop/QCDB/'
dat_names = ['snap_201801.txt', 'snap_201802-201803.txt', 'snap_201804-201805.txt']
my_dict = dict()

tot = pd.DataFrame()
for name in dat_names:
    my_dict[name] = pd.read_csv(dat_path+name, delim_whitespace=True, header=None)
    tot = tot.append(my_dict[name], ignore_index=True)
#%%

tot[2] = tot.index+1

#%%
## To replace with formatted printing to a file
tot.to_csv(dat_path+'db.arrival', index=False)

#%%

#test = pd.read_csv(dat_path+'db.arrival')