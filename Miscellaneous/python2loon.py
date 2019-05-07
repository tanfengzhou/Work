#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 11:23:38 2017

@author: visser
"""

import pandas as pd
import numpy as np
#import geopandas as gpd
from shapely.geometry import Point
import os
import obspy

#%%

'''
This script is designed to duplicate the results of db2loon, but not the format.
The final results includes the same columns as db2loon with a corrected magnitude.

The input files needed are the database path, and a table of moment magnitudes to
replace the ML values within the antelope database.

Six tables are needed from the Antelope db. 
These include: assoc, origin, origerr, stamag, arrival, netmag

origin, origerr, and netmag can be merged.
assoc and arrival can be merged.
stamag and assoc should be joined so that stamag can be sorted by time.

columns from the origin table that are needed are:
time,lat,lon,depth,ndef,ml.

columns from the origerr table that are needed are:
sdobs, stime, strike, sdepth, smajax, sminax, sdepth.

columns from the assoc table that are needed are:
timeres, phase, timedef, timeres, wgt, delta, esaz.

columns from the stamag table that are needed are:
sta, magtype, magnitude.

columns from the netmag table that are needed are:
magtype, magnitude, uncertainty, nsta.

!!! NEED TO HAVE FIXED ARRIVAL MAGNITUDES AND EVENT MAGNITUDES CALCULATED
BY THE MAGNITUDE CORRECTION SCRIPT !!!
    
'''

db_path = r'C:\Users\User\Desktop\Work Stuff\Scripts\data\pf299\\'
db_name = r'pf299'
auth = 'ISR  '
StDly = '0.00'
blank = ' '

canada_shp_file = gpd.read_file(r'C:\Users\User\Desktop\Work Stuff\Scripts\data\metadata\canada_names.csv.shp')
crs = {'init': 'epsg:4326'}

mw_events = pd.read_csv(r'C:\Users\User\Desktop\Work Stuff\Scripts\data\metadata\MW_catalogue.csv')
mw_events['date']=pd.to_datetime(mw_events['Y-M-D']+' '+mw_events['H:M:S'], format='%Y-%m-%d %H:%M:%S').dt.strftime('%Y-%m-%d %H:%M:%S')

#%%
site = pd.read_csv(r'C:\Users\User\Desktop\Work Stuff\Scripts\data\pf299\site.csv', header=None, usecols=[0,3,4], names=['sta','lat','lon'])
site = site.set_index('sta')

stamag = pd.read_csv(r'C:\Users\User\Desktop\Work Stuff\Scripts\data\pf299\pf299.stamag', usecols = [0,1,2,3,4,6,7,8], names=['magid','sta','arid','orid','evid','magtype','magnitude','ct'], delim_whitespace=True, header=None)

origin = pd.read_csv(db_path+db_name+'.origin', delim_whitespace=True, header=None, 
                     usecols = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23], 
                     names=['lat','lon','depth','time','orid','evid','jdate','nass','ndef','ndp','grn','srn','etype','review','depdp','dtype','mb','mbid','ms','msid','ml','mlid','algorithm','auth','comm','iddate'])
arrival = pd.read_csv(db_path+db_name+'.arrival', delim_whitespace=True, header=None, 
                     usecols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23],
                     names=['time','arid','jdate','stassid','chanid','chan','iphase','stype','deltim','azimuth','delaz','slow','delslo','ema','rect','amp','per','logat','clip','fm','snr','qual','auth'])
assoc = pd.read_csv(db_path+db_name+'.assoc', delim_whitespace=True, header=None, 
                     usecols = [0,1,2,3,4,5,6,7,8,9,10,15],
                     names=['arid','orid','sta','phase','belief','delta','seaz','esaz','timeres','timedef','azres','wgt'])
origerr = pd.read_csv(db_path+db_name+'.origerr', delim_whitespace=True, header=None, 
                     usecols = [0,11,12,13,14,15,16,], 
                     names=['orid','sdobs','smajax','sminax','strike','sdepth','stime'])

origerr = origerr[origerr.smajax<=10]

netmag = pd.read_csv(db_path+db_name+'.netmag', usecols = [2,3,4,5,6,7], names=['orid','evid','magtype','nsta','magnitude','ct'], delim_whitespace=True, header=None)

origerr_netmag = pd.merge(origerr, netmag, how='inner', on='orid')
origin_origerr_netmag = pd.merge(origerr_netmag,origin, how='inner', on='orid')

origin_origerr_netmag = origin_origerr_netmag[((origin_origerr_netmag.ml!=-999) & (origin_origerr_netmag.sdobs<1) & (origin_origerr_netmag.lat<61) & (origin_origerr_netmag.lat>52) & (origin_origerr_netmag.lon<-115) & (origin_origerr_netmag.lon>-126))]

#%%

#'DataFrame' object has no attribute 'mag'!!!

fixed_stamag = pd.read_csv('../../magnitude_correction/output_files/'+db_name+'_fixed_station_mags.csv', usecols = [3,8])
fixed_evmag = pd.read_csv('../../magnitude_correction/output_files/'+db_name+'_fixed_event_mags.csv', usecols = [0,1,2,3],index_col=[0])

#** FOLLOWING LINES CALCULATES EVENT MAGNITUDES FOR EACH EARTHQUAKES
stamag_unique = stamag.magid.unique()
event_mag = []

for magid in stamag_unique:
    var = np.median(stamag[stamag.magid==magid].magnitude).round(2)
    event_mag.append(np.around(var, decimals=2))
    
    stamag.loc[stamag.magid==magid, 'evmag'] = var
    
stamag['magdif'] = stamag.magnitude - stamag.evmag
#**

#%%

# =============================================================================
# CREATE A DYNAMIC CORRECTION FACTOR TABLE FROM ALL EARTHQUAKES IN THE DATABASE
# =============================================================================

cf_df = pd.DataFrame(columns=['station no', 'station name', 'correction factor (mean)', 'correction factor (median)', 'No. of events'])

cf_df['station name']=stamag.sta.unique()
cf_df = cf_df.sort_values('station name', ascending=True)
cf_df.reset_index(inplace=True,drop=True)
station_number=np.arange(1,len(cf_df)+1)
cf_df['station no']=station_number

def correction_factor_creator(sta_mag,cf_df):
    
    for station in cf_df['station name']:
        
        cf_number   = len(sta_mag[sta_mag['sta']==station])
        cf_median   = np.median(sta_mag[sta_mag['sta']==station].magdif)
        cf_mean     = np.mean(sta_mag[sta_mag['sta']==station].magdif)

        cf_df.loc[cf_df['station name']==station, 'No. of events'] = int(cf_number)
        cf_df.loc[cf_df['station name']==station, 'correction factor (median)'] = float(cf_median)
        cf_df.loc[cf_df['station name']==station, 'correction factor (mean)'] = float(cf_mean)
    
    return cf_df

cf_df = correction_factor_creator(stamag,cf_df)

cf_df.plot(x='No. of events', y='correction factor (mean)', marker='o', linewidth=0, title='Mean Correction Factors for Stations')

#%%

#Calculates a correction factor with provisions for hypocentral distance. Will probably loop over events

def mag_correction(input_ml, sta_lat, sta_lon, ev_lat, ev_lon, depth, s):
    temp_tuple = obspy.geodetics.base.gps2dist_azimuth(sta_lat, sta_lon, ev_lat, ev_lon)
    epi_dist = (temp_tuple[0])/1000
    d_hypo = np.sqrt((epi_dist**2) + depth**2)
    if (d_hypo<=85):
        output_ml = [(0.7974*np.log10(d_hypo/100)) + (0.0016*(d_hypo - 100)) + (-2.28566774*np.log10(epi_dist)+1.53265) + 3 + s + input_ml]
    else:
        output_ml = [(-0.1385*np.log10(d_hypo/100)) + (0.0016*(d_hypo - 100)) + (-2.28566774*np.log10(epi_dist)+1.53265) + 3 + s + input_ml]
        
    return output_ml

#%%
    
fixed_mag = []

for i in range(1, 2):
    temp_sta = stamag.sta[i]
    temp_orid = stamag.orid[i]
    temp_origin = origin[origin.orid==temp_orid]
    temp_cf_df = cf_df[cf_df['station name']==temp_sta]
    
    input_ml = stamag.magnitude[i]
    sta_lat = site.lat[temp_sta]
    sta_lon = site.lon[temp_sta]
    ev_lat = temp_origin.lat[temp_orid-1]
    ev_lon = temp_origin.lon[temp_orid-1]
    depth = temp_origin.depth[temp_orid-1]
    s = temp_cf_df['correction factor (mean)'][temp_cf_df.index[0]]
    
    output_ml = mag_correction(input_ml, sta_lat, sta_lon, ev_lat, ev_lon, depth, s)
    fixed_mag.append(output_ml)

#%%
#Needs a table of fixed magnitudes

stamag = pd.merge(stamag, fixed_stamag, how='inner', on='arid')
origin_origerr_netmag = pd.merge(origin_origerr_netmag, fixed_evmag, how='inner', on='orid')
del fixed_stamag, fixed_evmag

#==============================================================================
# Find number of stations used in earthquake locations for each ORID
#==============================================================================

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
    


origin_origerr_netmag = pd.merge(origin_origerr_netmag,number_of_stations, how='left', on='orid')

section1_cols = ['C','TF','orid','YearMoDy', 'HrMn', 'Secnd', 'Latitude', 'Longitude', 'Depth', '#St-', '#Ph','MT', 'Maggerss','Agncy']
#origin_origerr_netmag['C']='S'
origin_origerr_netmag['TF']='  '
origin_origerr_netmag['Agncy']=auth
origin_origerr_netmag['VM']='01'
origin_origerr_netmag['Weight']='WT ON '
origin_origerr_netmag['L']='N'

section1_1_cols = ['C','VM','L','Weight','RMS-','TErr-','LatErr--','LonErr--', 'DErr-','MajE','MinE','VerE','AzHor','Agncy']
#C  VM  L  Weight  RMS  TErr-  LatErr--  LonErr---  DErr-  MajE  MinE  VerE  AzHor Agncy
#%%

geometry = [Point(xy) for xy in zip(origin_origerr_netmag.lon, origin_origerr_netmag.lat)]
gdf = gpd.GeoDataFrame(origin_origerr_netmag, crs = crs, geometry=geometry)

alberta = canada_shp_file[canada_shp_file.NAME=='ALBERTA, CANADA']
bc = canada_shp_file[canada_shp_file.NAME=='BRITISH COLUMBIA, CANADA']
yukon = canada_shp_file[canada_shp_file.NAME=='YUKON TERRITORY, CANADA']
nwt = canada_shp_file[canada_shp_file.NAME=='NORTHWEST TERRITORIES, CANADA']
#%%
ab_df = gpd.tools.sjoin(gdf,alberta, how='inner')
bc_df = gpd.tools.sjoin(gdf,bc, how='inner')
ny_df = gpd.tools.sjoin(gdf,yukon, how='inner')
nwt_df = gpd.tools.sjoin(gdf,nwt, how='inner')

df = pd.concat([ab_df, bc_df, ny_df, nwt_df])

#%%
df['Date'] = pd.to_datetime(df.time, unit='s').dt.strftime('%Y-%m-%d %H:%M:%S')
df['YrMnDy'] = pd.to_datetime(df.time, unit='s').dt.strftime('%Y%m%d')
df['hrmn'] = pd.to_datetime(df.time, unit='s').dt.strftime('%H%M')
df['s'] = np.round(pd.to_datetime(df.time, unit='s').dt.strftime('%S.%f').astype(float),2).apply(lambda x: '{0:0.2f}'.format(x))
df['s'] = df['s'].str.zfill(5)
df.nsta = df.nsta.astype(str)

df = df.sort_values(by='Date', ascending=True)
df.reset_index(drop=True,inplace=True)

for i in range(0, len(df)):
    strike = df.strike.iloc[i]
    majaxerr = df.smajax.iloc[i]
    
    df.set_value(i, 'LatErr', abs( majaxerr / 2.0*np.cos( strike * 2.0 * 3.1416 / 360.0 ) ))
    df.set_value(i, 'LonErr', abs( majaxerr / 2.0*np.sin( strike * 2.0 * 3.1416 / 360.0 ) ))
    
#df['LatErr'] = 
#df['LonErr'] = abs( df.smajax / 2.0*np.sin( df.strike * 2.0 * 3.1416 / 360.0 ) )


for i in range(0,len(mw_events)):
    matching_eq_index_val = df[df.Date==mw_events.date.loc[i]].index.get_values()
    
    df.set_value(matching_eq_index_val, 'lat', mw_events['Lat. (N)'].loc[i])
    df.set_value(matching_eq_index_val, 'lon', mw_events['Lon. (E)'].loc[i])
    df.set_value(matching_eq_index_val, 'mag', mw_events.mag.loc[i])
    df.set_value(matching_eq_index_val, 'fixed_mag', mw_events.mag.loc[i])
    df.set_value(matching_eq_index_val, 'magtype', 'MW')
#    df.set_value(matching_eq_index_val, 'ct', 0)
#    df.set_value(matching_eq_index_val, 'nsta', '')
    

#%%

solution = pd.DataFrame(columns=section1_cols)
error = pd.DataFrame(columns=section1_1_cols)
magnitude = pd.DataFrame()
#%%

def input_values(input_table,output_table,column_in,column_out):
    output_table[column_out] = input_table[column_in]
    return output_table

section1_cols = ['C','TF','orid','YearMoDy', 'HrMn', 'Secnd', 'Latitude', 'Longitude', 'Depth', '#St-', '#Ph','MT', 'Maggerss','Agncy']

solution = input_values(df,solution,'orid','orid')
solution = input_values(df,solution,'YrMnDy','YearMoDy')
solution = input_values(df,solution,'hrmn','HrMn')
solution = input_values(df,solution,'s','Secnd')
solution = input_values(df,solution,'lat','Latitude')
solution = input_values(df,solution,'lon','Longitude')
solution = input_values(df,solution,'depth','Depth')
solution = input_values(df,solution,'sta_num','#St-')
solution = input_values(df,solution,'ndef','#Ph')
solution = input_values(df,solution,'magtype','MT')
solution = input_values(df,solution,'fixed_mag','Maggerss')
solution['TF'] = '  '
solution['C'] = 'S'
solution['Agncy'] = auth
solution['MT'] = solution['MT'].str.upper()



error = input_values(df,error,'orid','orid')
error = input_values(df,error,'VM','VM')
#section1 = input_values(df,section1,'TF','TF')
error = input_values(df,error,'Weight','Weight')
error = input_values(df,error,'sdobs','RMS-')
error = input_values(df,error,'stime','TErr-')
error = input_values(df,error,'LatErr','LatErr--')
error = input_values(df,error,'LonErr','LonErr--')
error = input_values(df,error,'sdepth','DErr-')
error = input_values(df,error,'smajax','MajE')
error = input_values(df,error,'sminax','MinE')
error = input_values(df,error,'sdepth','VerE')
error = input_values(df,error,'strike','AzHor')
error = input_values(df,error,'Agncy','Agncy')
error['C'] = 'E'
error['L'] = 'N'
error['VerE']=error['VerE'].astype(float)
error['AzHor']=error['AzHor'].astype(float)

magnitude['C'] = 'M'
magnitude = input_values(df,magnitude,'orid','#orid')
magnitude = input_values(df,magnitude,'magtype','magtype')
magnitude = input_values(df,magnitude,'magnitude','ML')
magnitude = input_values(df,magnitude,'fixed_mag','CML')
magnitude['MW']=''
magnitude = input_values(df,magnitude,'ct','ct')
magnitude = input_values(df,magnitude,'nsta','#StM')
                         

mw_magnitude = pd.DataFrame()
#mw_magnitude['MW']=''
#mw_magnitude['orid']=''

for i in range(0,len(mw_events)):
    matching_eq_index_val = df[df.Date==mw_events.date.loc[i]].index.get_values()
    
    mw_magnitude.set_value(i, 'MW', mw_events['mag'].loc[i])
    mw_magnitude.set_value(i, 'orid', df['orid'].loc[matching_eq_index_val[0]])
        
mw_magnitude['magtype'] = 'MW'

#%%
def number_of_decimals(input_table, column_name, number_of_decimals):

    if number_of_decimals ==0:
        input_table[column_name] = np.round(input_table[column_name], 0).apply(lambda x: '{0:0.0f}'.format(x))
    
    if number_of_decimals ==1:
        input_table[column_name] = np.round(input_table[column_name], 1).apply(lambda x: '{0:0.1f}'.format(x))
    
    if number_of_decimals ==2:
        input_table[column_name] = np.round(input_table[column_name], 2).apply(lambda x: '{0:0.2f}'.format(x))

    if number_of_decimals ==3:
        input_table[column_name] = np.round(input_table[column_name], 3).apply(lambda x: '{0:0.3f}'.format(x))

    if number_of_decimals ==4:
        input_table[column_name] = np.round(input_table[column_name], 4).apply(lambda x: '{0:0.4f}'.format(x))

    return input_table
    
solution = number_of_decimals(solution, 'Latitude', 4)
solution = number_of_decimals(solution, 'Longitude', 4)
solution = number_of_decimals(solution, 'Depth', 2)
solution = number_of_decimals(solution, 'Maggerss', 2)
error = number_of_decimals(error, 'LatErr--', 3)
error = number_of_decimals(error, 'LonErr--', 3)
error = number_of_decimals(error, 'TErr-', 2)
error = number_of_decimals(error, 'RMS-', 2)
error = number_of_decimals(error, 'DErr-', 2)
error = number_of_decimals(error, 'MajE', 2)
error = number_of_decimals(error, 'MinE', 2)
error = number_of_decimals(error, 'VerE', 2)
error = number_of_decimals(error, 'AzHor', 1)
error['AzHor'] = error['AzHor'].astype(str)
error['AzHor'] = error['AzHor'].str.zfill(5)
#%%
def align_strings(input_table,input_column):
    col_length = len(input_column)
    max_len = input_table[input_column].astype(str).map(len).max()
    
    difference = col_length-max_len
    
    spaces=difference*' '
    
    try:
        input_table[input_column].astype(float)
        
        input_table[input_column] = spaces + input_table[input_column].astype(str)
        
        print(str(len(spaces))+' added')
        
    except:
        print('column must be a string, so aligning left')
        
        try:
            input_table[input_column] = input_table[input_column].astype(str)+spaces
            
        except:
            pass
    
    return input_table
    print(input_table[input_column])
    
solution = align_strings(solution,'Latitude')
solution = align_strings(solution,'Longitude')
#solution = align_strings(solution,'Depth')
       #%%                          
error['LatErr--'] = error['LatErr--'].astype(str)+'km'
error['LonErr--'] = error['LonErr--'].astype(str)+'km'

def mixed_len_column_align_right(input_table,input_column):
    col_length = len(input_column)
    min_length = input_table[input_column].astype(str).map(len).min()
    
    for i in range(min_length, col_length):
        
        difference = col_length-i
        spaces=difference*' '
        
        input_table[input_column].ix[input_table[input_column].astype(str).map(len) == i] = spaces+input_table[input_column].astype(str)
        
#        input_table[input_column].ix[]
#    difference = col_length-max_len
    
#    for 
#    print(mask, input_table[input_column])
    return input_table
solution = mixed_len_column_align_right(solution,'Depth')
solution = mixed_len_column_align_right(solution,'#St-')
solution = mixed_len_column_align_right(solution,'#Ph')
solution = mixed_len_column_align_right(solution,'Maggerss')
error = mixed_len_column_align_right(error,'TErr-')
error = mixed_len_column_align_right(error,'LatErr--')
error = mixed_len_column_align_right(error,'LonErr--')
error = mixed_len_column_align_right(error,'DErr-')
error = mixed_len_column_align_right(error,'AzHor')
magnitude = mixed_len_column_align_right(magnitude,'ct')
magnitude = mixed_len_column_align_right(magnitude,'#StM')
#%%

section2 = pd.DataFrame()

section2['orid'] = df.orid
section2['C'] = 'C'
section2['Location (English)'] = df.NAME
section2['Location (French)'] = df.NOM
section2['blank'] = '    '
section2['Agncy'] = auth

#%%

def mixed_len_column_align_left(input_table,input_column):
    col_length = input_table[input_column].astype(str).map(len).max()
    min_length = input_table[input_column].astype(str).map(len).min()
    
    for i in range(min_length, col_length):
        
        difference = col_length-i
        spaces=difference*' '
        
        input_table[input_column].ix[input_table[input_column].astype(str).map(len) == i] = input_table[input_column].astype(str)+ spaces
        
    return input_table




section2 = mixed_len_column_align_left(section2,'Location (English)')
section2 = mixed_len_column_align_left(section2,'Location (French)')

#solution['Depth']=solution['Depth'].astype(float)
#%%
magnitude = number_of_decimals(magnitude, 'ct', 2)
magnitude = number_of_decimals(magnitude, 'ML', 2)
magnitude = number_of_decimals(magnitude, 'CML', 2)
magnitude['C'] = 'M'
magnitude['magtype'] = magnitude.magtype.str.upper()

error['DErr-'].ix[df['dtype']=='g'] = 'F    '


#%%

section3_cols = ['C','Statn', 'IC', 'nHHMM', 'SSSSS', 'TCorr','Q-Phase', 'IUW', 'TTres', 'LocW','StDly', 'EDistnc','Azm','Agncy']
section3 = pd.DataFrame(columns=section3_cols)

assoc_arrival = pd.merge(arrival, assoc, how='inner', on='arid')

assoc_arrival['IC'] = assoc_arrival.chan.str[0:1]+assoc_arrival.chan.str[2:3]

assoc_arrival['Date'] = pd.to_datetime(assoc_arrival.time, unit='s').dt.strftime('%Y-%m-%d %H:%M:%S')
assoc_arrival['YrMnDy'] = pd.to_datetime(assoc_arrival.time, unit='s').dt.strftime('%Y%m%d')
assoc_arrival['hrmn'] = pd.to_datetime(assoc_arrival.time, unit='s').dt.strftime('%H%M')
assoc_arrival['s'] = np.round(pd.to_datetime(assoc_arrival.time, unit='s').dt.strftime('%S.%f').astype(float),2).apply(lambda x: '{0:0.2f}'.format(x))
assoc_arrival['s'] = assoc_arrival['s'].str.zfill(5)

assoc_arrival = assoc_arrival[assoc_arrival.orid.isin(df.orid)]

section3['TCorr']='0.00'
section3['StDly']='0.00'
#assoc_arrival['station'] = assoc_arrival['sta_x']



assoc_arrival['timedef'].ix[assoc_arrival['timedef']!='d'] = 'x0'
assoc_arrival['timedef'].ix[assoc_arrival['timedef']=='d'] = '0'

assoc_arrival['edistance'] = assoc_arrival.delta*111.195

assoc_arrival = assoc_arrival.sort_values(by='edistance', ascending=True)
assoc_arrival.reset_index(drop=True,inplace=True)

section3 = input_values(assoc_arrival,section3,'orid','orid')
section3 = input_values(assoc_arrival,section3,'sta','Statn')
section3 = input_values(assoc_arrival,section3,'IC','IC')
section3 = input_values(assoc_arrival,section3,'hrmn','nHHMM')
section3 = input_values(assoc_arrival,section3,'s','SSSSS')
section3 = input_values(assoc_arrival,section3,'phase','Q-Phase')
section3 = input_values(assoc_arrival,section3,'timedef','IUW')
section3 = input_values(assoc_arrival,section3,'timeres','TTres')
section3 = input_values(assoc_arrival,section3,'wgt','LocW')
section3 = input_values(assoc_arrival,section3,'edistance','EDistnc')
section3 = input_values(assoc_arrival,section3,'esaz','Azm')
section3['Agncy'] = auth

section3['Statn'].ix[section3['Statn'].str.len()==4] = section3['Statn']+' '
section3['Statn'].ix[section3['Statn'].str.len()==3] = section3['Statn']+'  '

section3 = number_of_decimals(section3, 'TTres', 2)
section3 = number_of_decimals(section3, 'LocW', 2)
section3 = number_of_decimals(section3, 'EDistnc', 1)
section3 = number_of_decimals(section3, 'Azm', 0)
section3['Azm'] = section3['Azm'].str.zfill(3)
section3['C']= blank
section3 = mixed_len_column_align_left(section3,'Q-Phase')
section3['Q-Phase']=' '+section3['Q-Phase']+'     '
section3['Agncy']='  '+section3['Agncy']
section3['TCorr']=' 0.00'
section3['StDly']=' 0.00'
#%%

#==============================================================================
# Section 4: Station magnitudes
#==============================================================================

assoc_arrival_stamag = pd.merge(assoc_arrival, stamag, how='inner', on=['orid','sta'])

assoc_arrival_stamag['edistance'] = assoc_arrival_stamag.delta*111.195

assoc_arrival_stamag = assoc_arrival_stamag.drop_duplicates(subset=['orid','sta'], keep='last')

assoc_arrival_stamag['Date'] = pd.to_datetime(assoc_arrival_stamag.time, unit='s').dt.strftime('%Y-%m-%d %H:%M:%S')
assoc_arrival_stamag['YrMnDy'] = pd.to_datetime(assoc_arrival_stamag.time, unit='s').dt.strftime('%Y%m%d')
assoc_arrival_stamag['hrmn'] = pd.to_datetime(assoc_arrival_stamag.time, unit='s').dt.strftime('%H%M')
assoc_arrival_stamag['s'] = np.round(pd.to_datetime(assoc_arrival_stamag.time, unit='s').dt.strftime('%S.%f').astype(float),2).apply(lambda x: '{0:0.2f}'.format(x))
assoc_arrival_stamag['s'] = assoc_arrival['s'].str.zfill(5)

assoc_arrival_stamag = assoc_arrival_stamag[assoc_arrival_stamag.orid.isin(df.orid)]

assoc_arrival_stamag = assoc_arrival_stamag.sort_values(by='edistance', ascending=True)
assoc_arrival_stamag.reset_index(drop=True,inplace=True)


#%%

section4_columns=['C','Statn','IC','nHHMM','SSSSS','TCorr','-Phase--','Period','-Amplitude--','orid','MT','Maggerss']

section4 = pd.DataFrame(columns=section4_columns)

section4 = input_values(assoc_arrival_stamag,section4,'orid','orid')
section4 = input_values(assoc_arrival_stamag,section4,'sta','Statn')
#section4 = input_values(assoc_arrival_stamag,section4,'magtype','MagType')
#section4 = input_values(assoc_arrival_stamag,section4,'magnitude','Mag')
section4 = input_values(assoc_arrival_stamag,section4,'fixed_mag','Maggerss')

section4['Statn'].ix[section4['Statn'].str.len()==4] = section4['Statn']+' '
section4['Statn'].ix[section4['Statn'].str.len()==3] = section4['Statn']+'  '

#section4['MagType'] = section4['MagType'].str.upper()

section4['Agncy'] = auth

section4 = number_of_decimals(section4, 'Maggerss', 2)
#section4 = number_of_decimals(section4, 'CMag', 2)
section4['C'] = 'A'
section4['MT']='   ML    '
section4.IC = '  '
section4.nHHMM = '     '
section4.SSSSS = '     '
section4.TCorr = '     '
section4['-Phase--'] = '        '
section4.Period = '     '
section4['-Amplitude--'] = '            '




#%%

def convert_type_to_int(input_table,column):
    input_table[column] = input_table[column].astype(int)
    return input_table
solution = convert_type_to_int(solution,'orid')
error = convert_type_to_int(error,'orid')
#section3 = convert_type_to_int(section3,'orid')
#section4 = convert_type_to_int(section4,'orid')

#%%


if os.path.exists('../output_files/pick_file.txt'):
    os.remove('../output_files/pick_file.txt')

with open('../output_files/pick_file.txt', 'a') as the_file:
    for line in solution.orid:
#        section1_line = solution.ilocsolution.head(0)]
#        section2_line = error[error.head(0)]
        section1_line = solution[solution.orid==line]
        section2_line = section2[section2.orid==line]
        section3_line = error[error.orid==line]        
        section5_line = section3[section3.orid==line]
        section4_line = section4[section4.orid==line]
#        
        section1_line = section1_line.drop('orid', axis=1)
        section2_line = section2_line.drop('orid', axis=1)
        section3_line = section3_line.drop('orid', axis=1)
        section4_line = section4_line.drop('orid', axis=1)
        section5_line = section5_line.drop('orid', axis=1)

        the_file.write('C  TF  YearMoDy  HrMn  Secnd  Latitude  Longitude  Depth  #St-  #Ph  -Magnitude--  Agncy'+'\n')
        line1_of_catalogue = section1_line.to_string(header=False,index=False, col_space=0)
        the_file.write(line1_of_catalogue+'\n')
        the_file.write('C  VM  L  Weight  RMS-  TErr-  LatErr--  LonErr--  DErr-  MajE  MinE  VerE  AzHor  Agncy'+'\n')
        line3_of_catalogue = section3_line.to_string(header=False,index=False, col_space=0)
        the_file.write(line3_of_catalogue+'\n')
        
        
        any_mw = mw_magnitude[mw_magnitude.orid==line]
        maggy = magnitude[magnitude['#orid']==line]
        if len(any_mw) != 0:
            the_file.write('M  *MW    '+str(any_mw.MW.iloc[0])+ '                                                                     ISR \n')
            the_file.write('M  ML     '+str(maggy.CML.iloc[0])+'  ('+str(maggy.ct.iloc[0])+')  '+ str(maggy['#StM'].iloc[0])+'                                                          ISR \n')
            the_file.write('M  ML     '+str(maggy.ML.iloc[0])+'  ('+str(maggy.ct.iloc[0])+')  '+ str(maggy['#StM'].iloc[0])+'                                                          ISR \n')
            
        else:         
            the_file.write('M  *ML    '+str(maggy.CML.iloc[0])+'  ('+str(maggy.ct.iloc[0])+')  '+ str(maggy['#StM'].iloc[0])+'                                                       ISR \n')
            the_file.write('M  ML     '+str(maggy.ML.iloc[0])+'  ('+str(maggy.ct.iloc[0])+')  '+ str(maggy['#StM'].iloc[0])+'                                                       ISR \n')

            
        the_file.write('C  E   '+section2_line['Location (English)'].iloc[0]+section2_line['blank'].iloc[0]+'                                           ISR \n')        
        the_file.write('C  F   '+section2_line['Location (French)'].iloc[0]+section2_line['blank'].iloc[0]+'                                          ISR \n')   
        
#        the_file.write('C  Statn  IC  nHHMM  SSSSS  TCorr  Q-Phase  IUW  TTres  LocW  StDly  EDistnc  Azm  Agncy'+'\n')

        line5_of_catalogue = section5_line.to_string(header=True,index=False)
        the_file.write(line5_of_catalogue+'\n')
        
        the_file.write('C  Statn  IC nHHMM  SSSSS  -TCorr--  -Phase---  -Period  -Amplitude  -Magnitude--  Agncy'+'\n')
        
        line4_of_catalogue = section4_line.to_string(header=False,index=False)
        the_file.write(line4_of_catalogue+'\n \n')
#

#

#        
##        line2_of_catalogue = curated_merged_arrival_catalog[curated_merged_arrival_catalog.orid==i]
##        line2_of_catalogue = line2_of_catalogue.drop('orid', axis=1)
##        line2_of_catalogue = line2_of_catalogue.to_string(index=False)
##        the_file.write(line2_of_catalogue+'\n'+'\n')

