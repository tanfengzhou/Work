#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 09:05:02 2019

@author: pgcseismolab
"""

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from collections import Counter
import seaborn as sns

#%%

site = pd.read_csv('~/Documents/chet/dbmaster/2018/nbc.site',usecols=[0,1,2,3,4], names=['sta','ondate','offdate','lat','lon'],  header=None, delim_whitespace=True)
snap1 = pd.read_csv('~/Documents/chet/dbs/201704-201707/snap_201704-201707.txt', header=None, delim_whitespace=True)
snap2 = pd.read_csv('~/Documents/chet/dbs/201801-201803/snap_201801-201803.txt', header=None, delim_whitespace=True)
snap3 = pd.read_csv('~/Documents/chet/dbs/201804-201811/snap_201804-201811.txt', header=None, delim_whitespace=True)

stas = snap1[0].append(snap2[0].append(snap3[0]))

ustas = list(set(stas))


#%%

on_datetime = []
off_datetime = []
for time in site.ondate:
    on_datetime.append(datetime.strptime(str(time), '%Y%j').date())
    
for time in site.offdate:
    off_datetime.append(datetime.strptime(str(time), '%Y%j').date())
    
site['on_datetime'] = on_datetime
site['off_datetime'] = off_datetime

#%%

month_dict = dict()
start = datetime(2017, 1, 1)
end = datetime(2018, 12, 31) 
days = pd.date_range(start=start, end=end, periods=365*2)

for day in days:
    my_day = day.date()
    tmp = site[(my_day < site.off_datetime) & (my_day > site.on_datetime) & (site.sta.isin(ustas))]
    month_dict[day] = len(tmp)

#%%
    
#plt.bar(list(month_dict.keys()), month_dict.values())
#plt.show()


#%%

day_dict = dict()
day_dict['stations'] = month_dict
a = pd.DataFrame.from_dict(day_dict)

#%%

#sns.barplot(data=a, x=a.index, y=a.stations, color='lightblue')

#%%

plt.plot(list(month_dict.keys()), month_dict.values())

#%%

active_days = pd.date_range(start=datetime(2017, 4, 1), end=datetime(2017, 7, 1))
mid = pd.date_range(start=datetime(2018, 1, 1), end=datetime(2018, 11, 23))
end = pd.date_range(start=datetime(2018, 12, 7), end=datetime(2018, 12, 31))
all_active = active_days.append(mid.append(end))

#%%

active_df = pd.DataFrame()
the_dict = dict()

for day in days:
    my_day = day.date()
    tmp = site[(my_day < site.off_datetime) & (my_day > site.on_datetime) & (site.sta.isin(ustas)) & (my_day in all_active)]
    #the_dict[day] = len(tmp)
    test = len(tmp)
    if test > 0:
        the_dict[day] = test
    else:
        the_dict[day] = np.nan

#plt.xticks(rotation=45)    
##plt.bar(list(the_dict.keys()), the_dict.values(), width=1, alpha=0.5, color='green')
#plt.show()
#plt.plot(list(the_dict.keys()), the_dict.values(), color='red')
#plt.show()

#%%

my_dict = dict()
my_dict['stations'] = the_dict
b = pd.DataFrame.from_dict(my_dict)
#%%

test1 = site[(site.off_datetime > datetime(2017, 7, 1).date()) & (site.on_datetime < datetime(2017, 4, 1).date()) & (site.sta.isin(ustas))]
test2 = site[(site.off_datetime > datetime(2018, 1, 1).date()) & (site.on_datetime < datetime(2018, 11, 23).date()) & (site.sta.isin(ustas))]
test3 = site[(site.off_datetime > datetime(2018, 12, 7).date()) & (site.on_datetime < datetime(2019, 1, 1).date()) & (site.sta.isin(ustas))]

#begin = b[(b.index < datetime(2017, 7, 1).date())]
#mid = b[(b.index < datetime(2018, 11, 23)) & (b.index > datetime(2018, 1, 1))]
#end = b[(b.index < datetime(2018, 9, 1)) & (b.index > datetime(2018, 12, 7))]

#plt.bar(begin.index, begin.stations)

#%%

the_list = []
for row in b.index:
    rowd = row.date()
    if rowd <= datetime(2017, 7, 1).date():
        the_list.append('b')
        
    elif ((rowd <= datetime(2018, 11, 23).date()) & (rowd >= datetime(2018, 1, 1).date())):
        the_list.append('m')
        
    elif ((rowd >= datetime(2018, 12, 7).date())):
        the_list.append('e')
        
    else:
        the_list.append('o')

b['categorical'] = the_list

#%%

#sns.barplot(data=b, x=b.index, y = b.stations, hue=b.categorical)


#%%

begin = b[b.categorical=='b']
middle = b[b.categorical=='m']
end = b[b.categorical=='e']

plt.xticks(rotation=45)
plt.bar(begin.index, begin.stations, width=1)
plt.bar(middle.index, middle.stations, width=1)
plt.bar(end.index, end.stations, width=1)
plt.xlabel('Time')
plt.ylabel('Number of stations')
#plt.legend()
plt.show()




















