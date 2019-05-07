#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 15:10:03 2018

@author: chet
"""
#%%

'''
Reads in an Antelope arrival table and outputs an Antelope readable arrival table
'''

import pandas as pd
input_file = r'C:\Users\User\Desktop\Work Stuff\data\cat.arrival'
output_file = r'C:\Users\User\Desktop\Work Stuff\data\cat_remade.arrival'

#%%

arr = pd.read_csv(input_file, delim_whitespace=True, header=None)
cat = arr.drop_duplicates(subset=[0,1,3,6,7])
cat.columns = ['sta','time','arid','jdate','stassid','chanid','chan','iphase','stype','deltim','azimuth','delaz','slow','delslo', 'ema', 'rect', 'amp', 'per', 'logat', 'clip', 'fm', 'snr', 'qual', 'auth', 'commid', 'Iddate']
cat = cat[['sta','chan','iphase','arid','time']]
cat = cat.reset_index(drop=True)
arid = list(cat.index+1)

cat.sta = cat.sta.astype(str)
cat.chan = cat.chan.astype(str)
cat.iphase = cat.iphase.astype(str)
cat.arid = cat.arid.astype(int)
cat.time = cat.time.astype(float)
#%%

f = open(output_file, 'w')

for n in range(len(cat.arid)):
    f.write('%-6s %17.5f %8d %8d %8d %8d %-8s %-8s %s %6.3f %7.2f %7.2f %7.2f %7.2f %7.2f %7.3f %10.1f %7.2f %7.2f %s %-2s %10d %s %-16s %7d %17.5f\n' % (cat.sta[n], cat.time[n], cat.arid[n], -1, -1, -1, cat.chan[n], cat.iphase[n], '-', -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -999.00, '-', '-', -1, '-', 'dbp:cnn', -1, -999999999.999))

f.close() 

#%%
