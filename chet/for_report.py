# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 10:17:23 2018

@author: User
"""

#%%

import pandas as pd
import b_value
import datetime as dt
import obspy


#%%

#Read in all tables as pandas dataframes

nedb_cat = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\nedb_cat.csv')

april1_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt1\var_april.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
april1_origerr = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt1\var_april.origerr', delim_whitespace=True, header=None, usecols = [0,11,12,13,14,15], names=['orid','sdobs','smajax','sminax','strike','sdepth'])

april2_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt2\var_april_after_12.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
april2_origerr = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\aprilpt2\var_april_after_12.origerr', delim_whitespace=True, header=None, usecols = [0,11,12,13,14,15], names=['orid','sdobs','smajax','sminax','strike','sdepth'])

mar_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\march\var_march.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
mar_origerr = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\march\var_march.origerr', delim_whitespace=True, header=None, usecols = [0,11,12,13,14,15], names=['orid','sdobs','smajax','sminax','strike','sdepth'])

feb_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\feb\feb99work.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
feb_origerr = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\feb\feb99work.origerr', delim_whitespace=True, header=None, usecols = [0,11,12,13,14,15], names=['orid','sdobs','smajax','sminax','strike','sdepth'])

jan_origin = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\jan\pf299.origin', header=None, delim_whitespace=True, usecols=[0,1,2,3,4,7,8,20], names=['lat','lon','depth','time','orid','nassoc','ndef','ml'])
jan_origerr = pd.read_csv(r'C:\Users\User\Desktop\Work\Scripts\data\jan\pf299.origerr', delim_whitespace=True, header=None, usecols = [0,11,12,13,14,15], names=['orid','sdobs','smajax','sminax','strike','sdepth'])

#%%
## Merge all the tables together

april1_origin_origerr = pd.merge(april1_origin, april1_origerr, on='orid')
april2_origin_origerr = pd.merge(april2_origin, april2_origerr, on='orid')
mar_origin_origerr = pd.merge(mar_origin, mar_origerr, on='orid')
feb_origin_origerr = pd.merge(feb_origin, feb_origerr, on='orid')
jan_origin_origerr = pd.merge(jan_origin, jan_origerr, on='orid')
cat = jan_origin_origerr.append(feb_origin_origerr.append(mar_origin_origerr.append(april1_origin_origerr.append(april2_origin_origerr))))

#%%
## Restrict data to study area

lat_max = 61
lat_min = 52
lon_max = -115
lon_min = -126
max_depth = 20
max_smajax = 10
max_sdobs = 1
min_ml = -999

cat = cat[(cat.lat<=lat_max)&(cat.lat>=lat_min)&(cat.lon<=lon_max)&(cat.lon>=lon_min)&(cat.depth<=max_depth)&(cat.smajax<=max_smajax)&(cat.sdobs<=max_sdobs)&(cat.ml>-999)]

#%%
##Reset the index

cat = cat.reset_index(drop=True)
cat.orid = cat.index + 1

#%%
## Summary statistics

cat_described = cat.describe()

#%%
##Create a column containing epoch time

nedb_cat['strepoch'] = nedb_cat.Date + ' ' + nedb_cat['Time(UT)']
nedb_cat['date'] = pd.to_datetime(nedb_cat.strepoch, format='%d/%m/%Y %H:%M:%S')
nedb_cat['epoch'] = (nedb_cat['date'] - dt.datetime(1970,1,1)).dt.total_seconds()

  
#%%
    
left = pd.DataFrame()
right = pd.DataFrame()
mergeid = 0
    
for i in range(len(cat)):
    for j in range(len(nedb_cat)):
        temp_dif = abs(cat.time[i]-nedb_cat.epoch[j])
        if(temp_dif<=8):
            temp_cat = cat.iloc[i]
            temp_cat['mergeid'] = mergeid
            temp_nedb = nedb_cat.iloc[j]
            temp_nedb['mergeid'] = mergeid
            left = left.append(temp_cat)
            right = right.append(temp_nedb)
            mergeid = mergeid + 1
            
to_compare = pd.merge(left, right, on='mergeid')


#%%

dist_list = []

for x in range(len(to_compare)):
    lat1 = to_compare.lat[x]
    lon1 = to_compare.lon[x]
    lat2 = to_compare.Lat[x]
    lon2 = to_compare.Long[x]
    
    temp_tuple = obspy.geodetics.base.gps2dist_azimuth(lat1, lon1, lat2, lon2)
    dist = temp_tuple[0]/1000
    dist_list.append(dist)    
    
to_compare['dist'] = dist_list

#%%

nedb_mags = nedb_cat.Mag.str.slice(start=0,stop=3).astype(float)

nedb_mc = b_value.find_best_fit(nedb_mags, 0, 3, 0.1)

#%%

my_mc = b_value.find_best_fit(cat.ml, 0, 3, 0.1)