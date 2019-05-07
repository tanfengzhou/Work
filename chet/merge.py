#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 15:46:56 2019

@author: pgcseismolab
"""

#%%

import pandas as pd


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
#tot.to_csv(dat_path+'db.arrival', index=False)

#%%

## Need to change variables for given file formats

output_file = dat_path+'db.arrival'
f = open(output_file, 'w')

for n in range(len(tot.index)):
    f.write('%-6s %17.5f %8d %8d %8d %8d %-8s %-8s %s %6.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %10.1f %7.2f %7.2f %s %-2s %10d %s %-16s %7d %17.5f\n' % (cat.sta[n], cat.time[n], cat.arid[n], -1, -1, -1, cat.chan[n], cat.iphase[n], '-', -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -999.00, '-', '-', -1, '-', 'dbp:cnn', -1, -999999999.999))

f.close() 
