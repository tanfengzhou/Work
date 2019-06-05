#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:42:52 2019

@author: pgcseismolab
"""
#%%

from obspy.clients.fdsn import Client
from obspy.core.utcdatetime import UTCDateTime
from obspy import read_inventory, read

import os
import time

#%%


#inventory = client.get_stations(station='AKSS', level='response', starttime = UTCDateTime("2016-11-13T11:59"), endtime=UTCDateTime("2016-11-13T13:01"))
#channel = inventory[0][0][0]
#resp = channel.response
##resp.plot()
#st_dict['AKSS'].remove_response(inventory=inventory, plot=True)

#%%

home = os.path.expanduser('~')
items = os.listdir(home + '/Documents/nz_antelope/wfs')
sta_names = []

for names in items:
    if names.endswith('.mseed'):
        sta_names.append(names.split('.')[0])
        
#%%

client = Client("GEONET")
for sta in sta_names:
    tmp_inv = client.get_stations(station=sta, level='response', starttime = UTCDateTime("2016-11-13T11:59"), endtime=UTCDateTime("2016-11-13T13:01"))
    tmp_inv.write(home + '/Documents/nz_antelope/station_xml/' + sta + '.xml', format='STATIONXML')
    time.sleep(1)
    
#%%

inv = read_inventory(home + '/Documents/nz_antelope/station_xml/AKSS.xml', format='STATIONXML')
st = read(home + '/Documents/nz_antelope/wfs/AKSS.mseed')
st.remove_response(inventory = inv)