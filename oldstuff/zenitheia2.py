# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:33:28 2022

Zenith angle calculation. 

VASIKA KAI POU EKANA THN METATROPH SE UTC+0 ARXIKA PALI SE TIPIKH WRA MOU
EVGLE TA APOTELESMATA. OPOTE KANW APO THN ARXH SE TOPIKH WRA TA PANTA KAI META
STO GHI TESTS ARXEIO KANW TO INDEX APO UTC+2 SE UTC+0??? 

FAILED CONCEPT. ZENITHEIA3 IS WHERE THE MONEY IS AT. 

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


# inserting datetimes as a DatetimeIndex 
data_file = 'Solar_1min_2021.txt'
df = pd.read_csv(data_file, usecols=[0], index_col=[0], header=None, parse_dates=True)


# convert greek standard time (GMT+2) to UTC 
grc_w = pytz.timezone('Etc/GMT-2')     # Greek Standard Timezone (UTC+2)      
utc=pytz.timezone('UCT')     # UTC (no offset)    
df.index = df.index.tz_localize(grc_w).tz_convert(utc) 


# make a Datetime column from the DateTimeIndex 
df['Datetime'] = df.index.to_series() 
# day count for the year
df['Day Count'] = df['Datetime'].dt.dayofyear 
# edw isws thelei allagh to 366 tou disektou etous, an den pairnei N=366
# o typos gia thn apoklish tou hliou? 


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


'''
Kai ena perasma se .txt arxeio menei??? 
Kalytera na enswmatwsv ton parapanw kwdika sto quality control programma ths 
ghi pou etoimazw hdh??
Tha dw ti tha kanw...
'''
