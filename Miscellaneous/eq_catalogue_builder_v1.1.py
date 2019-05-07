#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 10:39:36 2017

@author: visser
"""

import pandas as pd
import numpy as np
#from pyspark.sql import functions as F
import os.path

#%%

#db_path = input('Please enter the path to your Antelope database:  ')
#db_name = input('Please enter the name of the database:  ')

db_path = '/home/chet/Documents/to_run_dbevproc/jan_roughmag_complete/'
db_name = 'pf299'

#%%

###### WARNING!!! MUST CRUNCH ALL TABLES AND RUN DBFIXIDS FOR ALL IDS ON DB #####

# =============================================================================
# LAT LON BOUNDS
# =============================================================================

lat_min = 52 #degrees North
lat_max = 61 #degrees North
lon_min = -126 #degrees East
lon_max = -115 #degrees East

# =============================================================================
# MAX ALLOWABLE ERROR ELLIPSE MAJOR AXIS DISTANCE
# =============================================================================

hypoellipse_max=10 # in km

# =============================================================================
# LOAD ANTELOPE TABLES
# =============================================================================

# LOAD ORIGIN TABLE
origin = pd.read_csv(db_path+db_name+'.origin', delim_whitespace=True, header=None, 
                     usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23], 
                     names=['lat','lon','depth','time','orid','evid','jdate','nass','ndef','ndp','grn','srn','etype','review','depdp','dtype','mb','mbid','ms','msid','ml','mlid','algorithm','auth','comm','iddate'])

origin = origin[origin.time>0] # removes nulls since times cannot be negative, and antelope's null value is -999
origin = origin[origin.orid>0] # removes remaining nulls if present

origin = origin[((origin.lat>=lat_min) & # Removes events outside of lat/lon bounds
                 (origin.lat<=lat_max) & 
                 (origin.lon<=lon_max) & 
                 (origin.lon>=lon_min))]

# ORIGIN FORMATTING
origin.reset_index(drop=True,inplace=True) 
origin['time'] = np.round(origin['time'], 1)
origin['date'] = pd.to_datetime(origin['time'], unit='s')
origin = origin.drop('time', axis=1)
origin['sta_num']=''

# LOAD ARRIVAL TABLE
arrival = pd.read_csv(db_path+db_name+'.arrival', delim_whitespace=True, header=None, 
                     usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],
                     names=['sta','time','arid','jdate','stassid','chanid','chan','iphase','stype','deltim','azimuth','delaz','slow','delslo','ema','rect','amp','per','logat','clip','fm','snr','qual','auth'])

# ARRIVAL FORMATTING
arrival['time'] = np.round(arrival['time'], 1)
arrival['date'] = pd.to_datetime(arrival['time'], unit='s')
arrival = arrival.drop('time', axis=1)

# LOAD ASSOC TABLE
assoc = pd.read_csv(db_path+db_name+'.assoc', delim_whitespace=True, header=None, 
                     usecols = [0,1,2,9],
                     names=['arid','orid','sta','timedef'])

# LOAD ORIGERR TABLE
origerr = pd.read_csv(db_path+db_name+'.origerr', delim_whitespace=True, header=None, 
                     usecols = [0,12,13,14,15], 
                     names=['orid','smajax','sminax','strike','sdepth'])

origerr = origerr[(origerr.smajax<=hypoellipse_max)] # Removes events with error greater than the hypoellipse_max parameter

# LOAD STAMAG TABLE
sta_mags = pd.read_csv(db_path+db_name+'.stamag', usecols = [0,1,2,3,4,7,8], names=['magid','sta','arid','orid','evid','mag','uncertainty'], delim_whitespace=True, header=None)

#** FOLLOWING LINES CALCULATES EVENT MAGNITUDES FOR EACH EARTHQUAKES
stamag_unique = sta_mags.magid.unique()
event_mag = []

for magid in stamag_unique:
    var = np.median(sta_mags[sta_mags.magid==magid].mag).round(2)
    event_mag.append(np.around(var, decimals=2))
    
    sta_mags.loc[sta_mags.magid==magid, 'evmag'] = var
    
sta_mags['magdif'] = sta_mags.mag - sta_mags.evmag
#**

# MERGE ARRIVAL AND ASSOC TABLES USING INNER JOIN ON ARID AND REMOVE NULL ROWS
arrival_merged = pd.merge(arrival,assoc, how='inner', on='arid') 
arrival_merged = arrival_merged[arrival_merged.orid>0]

# MERGE ORIGIN AND ORIGERR TABLES ON ORID USING AN INNER JOIN ON ORID
origin_merged = pd.merge(origin,origerr, how='inner', on='orid')

# ONLY KEEP STATION MAGNITUDES THAT EXIST WITHIN THE ORIGIN MERGED TABLE AFTER BEING FILTERED
sta_mags = sta_mags[sta_mags.orid.isin(origin_merged.orid)]

#%%

# =============================================================================
# CREATE A DYNAMIC CORRECTION FACTOR TABLE FROM ALL EARTHQUAKES IN THE DATABASE
# =============================================================================

cf_df = pd.DataFrame(columns=['station no', 'station name', 'correction factor (mean)', 'correction factor (median)', 'No. of events'])

cf_df['station name']=sta_mags.sta.unique()
cf_df = cf_df.sort_values('station name', ascending=True)
cf_df.reset_index(inplace=True,drop=True)
station_number=np.arange(1,len(cf_df)+1)
cf_df['station no']=station_number

def correction_factor_creator(stamag,cf_df):
    
    for station in cf_df['station name']:
        
        cf_number   = len(stamag[stamag['sta']==station])
        cf_median   = np.median(stamag[stamag['sta']==station].magdif)
        cf_mean     = np.mean(stamag[stamag['sta']==station].magdif)

        cf_df.loc[cf_df['station name']==station, 'No. of events'] = int(cf_number)
        cf_df.loc[cf_df['station name']==station, 'correction factor (median)'] = float(cf_median)
        cf_df.loc[cf_df['station name']==station, 'correction factor (mean)'] = float(cf_mean)
    
    return cf_df

cf_df = correction_factor_creator(sta_mags,cf_df)

cf_df.plot(x='No. of events', y='correction factor (mean)', marker='o', linewidth=0, title='Mean Correction Factors for Stations')

#%%

correction_factors_low_number_of_events = cf_df[cf_df['No. of events']<=5]

#%%

stamag_new = pd.DataFrame()

for i in range(0, len(cf_df)):
    station = cf_df['station name'].loc[i]
    correction_factor = cf_df['correction factor (mean)'].loc[i]
    stamag_filtered = sta_mags[sta_mags['sta']==station]
    stamag_filtered['fixed_mag'] = stamag_filtered.mag-correction_factor
    stamag_new = stamag_new.append(stamag_filtered)
    
stamag_new.reset_index(drop=True, inplace=True)
stamag_new = stamag_new.sort_values(['orid'], ascending=True)

unique_orids = stamag_new.orid.unique()

orid_lst=[]
orid_mag_original=[]
orid_mag_fixed = []
bad_correction_factor_used=[]

for orid in unique_orids:
    
    stamag_filtered=stamag_new[stamag_new['orid']==orid]
    stamag_filtered.reset_index(drop=True, inplace=True)
    orid_mag_original.append(np.median(stamag_filtered['mag']))
    orid_mag_fixed.append(np.round(np.median(stamag_filtered['fixed_mag']),2))
    orid_lst.append(stamag_filtered['orid'].loc[0])
    
    yes_or_no = stamag_filtered[stamag_filtered.sta.isin(correction_factors_low_number_of_events['station name'])]
    
    if len(yes_or_no)!=0:
        bad_correction_factor_used.append('Yes')
    else:
        bad_correction_factor_used.append('No')

event_magnitude=pd.DataFrame()

event_magnitude['orid']=orid_lst
event_magnitude['mag']=orid_mag_original
event_magnitude['fixed_mag']=orid_mag_fixed
event_magnitude['bad_factor_used']=bad_correction_factor_used

#%%

cols = ['date','lat','lon','orid']

missing_event = origin[~origin.orid.isin(event_magnitude.orid)]
filtered_missing_events = missing_event[cols]
#filtered_missing_events.to_csv('../output_files/missing_events.csv', index=False)

#%%

origin_merged = origin_merged[origin_merged.orid.isin(orid_lst)]

#%%

origin_merged = pd.merge(origin_merged, event_magnitude, on='orid')

#%%

origin_merged.fixed_mag = np.round(origin_merged.fixed_mag,2)
origin_merged.depth = np.round(origin_merged.depth,1)
origin_merged.ml = np.round(origin_merged.ml,1)
origin_merged.fixed_mag = np.round(origin_merged.fixed_mag,1)

origin_merged.reset_index(drop=True,inplace=True)

#%%

origin_merged_1 = origin_merged[origin_merged.bad_factor_used!='Yes']
origin_merged_2 = origin_merged[origin_merged.bad_factor_used=='Yes']
origin_merged_2['fixed_mag'] = origin_merged_2['fixed_mag'].astype(str)+'&'
origin_merged = pd.concat([origin_merged_1,origin_merged_2])
origin_merged = origin_merged[origin_merged['fixed_mag'].notnull()]

origin_merged.smajax = np.round(origin_merged.smajax,1)
origin_merged.sminax = np.round(origin_merged.sminax,1)
origin_merged.sdepth = np.round(origin_merged.sdepth,1)
origin_merged.strike = np.round(origin_merged.strike,1)
origin_merged.loc[origin_merged.dtype == 'g', 'sdepth'] = 'F'

mw_events = pd.read_csv('../../python2loon/input_files/MW_catalogue.csv')
mw_events['date']=pd.to_datetime(mw_events['Y-M-D']+' '+mw_events['H:M:S'], format='%Y-%m-%d %H:%M:%S').dt.strftime('%Y-%m-%d %H:%M:%S')

for i in range(0,len(mw_events)):
    matching_eq_index_val = origin_merged[origin_merged.date.dt.strftime('%Y-%m-%d %H:%M:%S')==mw_events.date.loc[i]].index.get_values()
    
    origin_merged.set_value(matching_eq_index_val, 'lat', mw_events['Lat. (N)'].loc[i])
    origin_merged.set_value(matching_eq_index_val, 'lon', mw_events['Lon. (E)'].loc[i])
#    origin_merged.set_value(matching_eq_index_val, 'mag', mw_events.mag.loc[i])
    origin_merged.set_value(matching_eq_index_val, 'MW', mw_events.mag.loc[i])
#    origin_merged.set_value(matching_eq_index_val, 'magtype', 'MW')
#    origin_merged.set_value(matching_eq_index_val, 'uncertainty', 0)
#    origin_merged.set_value(matching_eq_index_val, 'nsta', '')

#%%

origin_merged['date'] = pd.to_datetime(origin_merged['date'], infer_datetime_format=True)
origin_merged = origin_merged.sort_values(['date'], ascending=True)
origin_merged.reset_index(drop=True, inplace=True)


time_defining = assoc[(assoc.timedef=='d') & ((assoc.orid>0))]
number_of_stations = pd.DataFrame()
number_of_stations['orid'] = ''
number_of_stations['sta_num'] = ''
unique_orid = time_defining.orid.unique()

for i in range(0,len(unique_orid)):
    orid = unique_orid[i]
    sta_num = len(time_defining[time_defining.orid==orid].sta.unique())
    number_of_stations.set_value(i, 'orid', orid)
    number_of_stations.set_value(i, 'sta_num', sta_num)

#%%

variance_mag = []
variance_cmag = []

stamag_new_newer = stamag_new[stamag_new.orid.isin(origin_merged.orid)]

for line in stamag_new_newer.orid.unique():
    variance_mag.append(np.var(stamag_new_newer[stamag_new_newer.orid==line].mag))
    variance_cmag.append(np.var(stamag_new_newer[stamag_new_newer.orid==line].fixed_mag))
    
mag_var = np.mean(variance_mag)
cmag_var = np.mean(variance_cmag)

print('Magnitude variance was originally '+str(np.round(mag_var, 2))+' and was corrected to '+str(np.round(cmag_var ,2)))

#%%

eq_catalog = pd.DataFrame()
eq_catalog['Y-M-D'] = origin_merged.date.dt.strftime('%Y-%m-%d')
eq_catalog['H:M:S'] = origin_merged.date.dt.strftime('%H:%M:%S.%f').apply(lambda x: x[:-5])
eq_catalog['Lat. (N)'] = origin_merged.lat
eq_catalog['Lon. (E)'] = origin_merged.lon
eq_catalog['ML'] = origin_merged.ml
eq_catalog['CML'] = origin_merged.fixed_mag.astype(str)
eq_catalog['CML'] = eq_catalog['CML'].replace('&','*', regex=True)
eq_catalog['MW'] = origin_merged.MW
eq_catalog['Depth'] = origin_merged.depth
#eq_catalog['#Sta']=origin_merged.sta_num
#eq_catalog['#Phases']=origin_merged.ndef
eq_catalog['MajorAxis Err']=origin_merged.smajax
eq_catalog['MinorAxis Err']=origin_merged.sminax
eq_catalog['Azimuth']=origin_merged.strike
eq_catalog['Depth Err']=origin_merged.sdepth

#%%

event_catalog = pd.DataFrame()
event_catalog['EQ Num'] = origin_merged.index.values+1
event_catalog['Y-M-D'] = origin_merged.date.dt.strftime('%Y-%m-%d')
event_catalog['H:M:S'] = origin_merged.date.dt.strftime('%H:%M:%S.%f').apply(lambda x: x[:-5])
event_catalog['Lat. (N)'] = origin_merged.lat
event_catalog['Lon. (E)'] = origin_merged.lon
event_catalog['ML'] = origin_merged.ml.astype(str)+'Ml'
event_catalog['CML'] = origin_merged.fixed_mag.astype(str)+'Ml'
event_catalog['Depth'] = origin_merged.depth.astype(str)+'km'
event_catalog['orid']=origin_merged.orid
event_catalog=event_catalog.set_index('orid')
event_catalog['CML'] = event_catalog['CML'].replace('&Ml','Ml*', regex=True)

#%%

#arrival_catalog = pd.DataFrame()
#arrival_catalog['orid']=arrival_merged.orid
#arrival_catalog['date']=arrival_merged.date
#arrival_catalog['sta']=arrival_merged.sta
#arrival_catalog['arid']=arrival_merged.arid
#arrival_catalog['chan']=arrival_merged.chan
#arrival_catalog['iphase']=arrival_merged.iphase

arrival_catalog = pd.merge(arrival, stamag_new, how='inner', on='arid')
arrival_catalog = arrival_catalog[arrival_catalog.orid>0]
arrival_catalog2 = pd.merge(arrival, assoc, how='inner', on='arid')
arrival_catalog2 = arrival_catalog2[arrival_catalog2.orid>0]

curated_arrival_catalog = pd.DataFrame()

curated_arrival_catalog['sta'] = arrival_catalog.sta_x
curated_arrival_catalog['chan'] = arrival_catalog.chan
curated_arrival_catalog['phase'] = arrival_catalog.iphase
curated_arrival_catalog['time'] = arrival_catalog.date.dt.strftime('%H:%M:%S.%f')
curated_arrival_catalog['snr'] = np.round(arrival_catalog.snr.astype(float),2).apply(lambda x: '{0:0.2f}'.format(x))
curated_arrival_catalog['ml'] = arrival_catalog.mag.apply(lambda x: '{0:0.2f}'.format(x))
curated_arrival_catalog['cml'] = np.round(arrival_catalog.fixed_mag,2).apply(lambda x: '{0:0.2f}'.format(x))
curated_arrival_catalog['orid'] = arrival_catalog.orid
curated_arrival_catalog['date'] = arrival_catalog.date

curated_arrival_catalog2 = pd.DataFrame()

curated_arrival_catalog2['sta'] = arrival_catalog2.sta_x
curated_arrival_catalog2['chan'] = arrival_catalog2.chan
curated_arrival_catalog2['phase'] = arrival_catalog2.iphase
curated_arrival_catalog2['time'] = arrival_catalog2.date.dt.strftime('%H:%M:%S.%f')
#curated_arrival_catalog2['snr'] = arrival_catalog2.snr
curated_arrival_catalog2['orid'] = arrival_catalog2.orid
curated_arrival_catalog2['date'] = arrival_catalog2.date

merged_arrival_catalog = pd.concat([curated_arrival_catalog2, curated_arrival_catalog])

merged_arrival_catalog = merged_arrival_catalog.sort_values('date', ascending=True)
merged_arrival_catalog.reset_index(drop=True,inplace=True)

merged_arrival_catalog = merged_arrival_catalog.drop('date', axis=1)

curated_merged_arrival_catalog = pd.DataFrame()

curated_merged_arrival_catalog['sta'] = merged_arrival_catalog.sta
curated_merged_arrival_catalog['chan'] = merged_arrival_catalog.chan
curated_merged_arrival_catalog['phase'] = merged_arrival_catalog.phase
curated_merged_arrival_catalog['time'] = merged_arrival_catalog.time.str[:-5]
curated_merged_arrival_catalog['ml'] = merged_arrival_catalog.ml
curated_merged_arrival_catalog['cml'] = merged_arrival_catalog.cml
curated_merged_arrival_catalog['snr'] = merged_arrival_catalog.snr
curated_merged_arrival_catalog['orid'] = merged_arrival_catalog.orid
curated_merged_arrival_catalog = curated_merged_arrival_catalog.fillna('')
#curated_merged_arrival_catalog['sta'] = merged_arrival_catalog.sta

#merged_arrival_catalog = merged_arrival_catalog.drop('date', axis=1)

#%%

orid_list= origin_merged.orid.tolist()

#%%

eq_catalog.to_csv('../output_files/eq_catalog_10kmV1.1.csv', index=False)

#%%

#sdobs_6 = len(origin_merged[origin_merged.sdobs<=0.6])
#sdobs_1 = len(origin_merged[origin_merged.sdobs<=1])
#
#sdobs_total = len(origin_merged)
#
#percent_6 = (sdobs_6/sdobs_total)*100
#percent_1 = (sdobs_1/sdobs_total)*100


#%%

if os.path.exists('/home/chet/Documents/constrained_dbs/sta_magV1.1.txt'):
    os.remove('/home/chet/Documents/constrained_dbs/sta_magV1.1.txt')

with open('/home/chet/Documents/constrained_dbs/sta_magV1.1.txt', 'a') as the_file:
    for i in orid_list:
        line_of_catalogue = event_catalog.loc[[i]]
        line_of_catalogue = line_of_catalogue.to_string(header=False,index=False)
        the_file.write(line_of_catalogue+'\n')
        line2_of_catalogue = curated_merged_arrival_catalog[curated_merged_arrival_catalog.orid==i]
        line2_of_catalogue = line2_of_catalogue.drop('orid', axis=1)
        line2_of_catalogue = line2_of_catalogue.to_string(index=False)
        the_file.write(line2_of_catalogue+'\n'+'\n')
