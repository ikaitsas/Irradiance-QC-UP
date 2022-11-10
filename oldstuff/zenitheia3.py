# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:33:28 2022

Zenith angle calculation. 
The Datetime column in in Greek standard time (UTC+2). It is converted to 
UTC+0 at the end, in the index.

NOMIZW BRHKA AYTO POU THELW NA KANW!! 

@author: yiann
"""

import pandas as pd
import numpy as np 
import pytz 


# latitude (lat) and longtitude (lon)
lat = 38.29136389     # (>0 north & <0 south) 
lon = 21.78858333     # (>0 east & <0 west)
# longitude of location's standard merindian (prime merindian at 0)
lon_m = 30.0     # Greece is +2 hours from 0 (15deg per standard merindian)     


# inserting datetimes
data_file = 'Solar_1min_2021.txt'
df = pd.read_csv(data_file, usecols=[0], header=None, parse_dates=True)
df.columns = ['Datetime'] 


# convert 'Datetime' to datetime object (alliws den doulevei tipota) 
df['Datetime'] = pd.to_datetime(df['Datetime']) 
# day count for the year
df['Day Count'] = df['Datetime'].dt.dayofyear 


# local time (LT) calculation  
# number of hours in a day
df['Hours'] = df['Datetime'].dt.hour
# number of minutes in an hour
df['Minutes'] = df['Datetime'].dt.minute 
# LT in decimal format   
df['LT'] = df['Hours'] + df['Minutes']/60  


# equation (correction) of time in decimal format (ET Decimal)
df['beta'] = (360/364) * ( df['Day Count'] - 81 )   
df['ET'] = (9.87 * np.sin(np.deg2rad(2*df['beta'])) - 7.53 * np.cos(np.deg2rad(df['beta'])) - 1.5 * np.sin(np.deg2rad(df['beta'])))/60           


# longtitude correction (LC) calculation, in decimal hours format 
df['LC'] = ((-4) * (lon_m - lon))/60


# true solar time (TST) calculation in decimal format
df['TST'] = df['LT'] + df['ET'] + df['LC'] 
# making negative values positive (alla allazei kai h stoixhsh se sxesh me thn mera?)
df['TST'] = np.where( df['TST'] < 0, df['TST'] + 24, df['TST'])


# hour angle calculation (h - in degrees) 
df['h'] = ( df['TST'] - 12 ) * 15    


# declination of the sun calculation (delta - in degrees)
df['delta'] = 23.45 * np.sin(np.deg2rad( (360/365) * (df['Day Count'] + 284) )) 


# cozZ calculation - where Z is the zenith angle
df['cosZ'] = np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(df['delta'])) + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(df['delta'])) * np.cos(np.deg2rad(df['h']))         


# zenith angle (Z - in degrees) calculation 
df['Z'] =  np.rad2deg(np.arccos(df['cosZ'])) 


# make the index a UTC+0 DatetimeIndex (original/local Datetime is UTC+2) 
df.index = df['Datetime'] 
df.index.names = ['Datetime UTC']   
grc_std = pytz.timezone('Etc/GMT-2')     # Greek Winter Timezone (+2h offset)      
utc=pytz.timezone('UCT')     # UTC (no offset)    
df.index = df.index.tz_localize(grc_std).tz_convert(utc) 


'''
Kai ena perasma se .txt arxeio menei??? 
Kalytera na enswmatwsv ton parapanw kwdika sto quality control programma ths 
ghi pou etoimazw hdh?? 
'''
