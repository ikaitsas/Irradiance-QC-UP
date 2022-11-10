# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

O kwdikas kanei elegxo twn ghi timwn
H ghi antistoixei sthn sthlh 6 apo tis 0-13 (14 sthles to arxeio) 

Mallon elysa to provlhma me ta duplicate indexes. Kai ta dyo df apoktoun ta 
idia indexes, apo thn idia phgh, opote kai mallon den tha yparxei provlhma
sthn stoixhsh twn ghi timwn kai twn timwn zenitheias gwnias.  

To provlhma me ta duplicate indexes mporei na to lysw kai me kapoion allon
tropo, an einai pio apodotikos (p.x. kana append - pros to paron to afhnw).
Etsi kai alliws kai ta duo auta indexes apo thn sthlh hmeromhnias kai wras
tou arxeiou tou aktinometrikou stathmou proekypsan (idia phgh), prin to 
dropna() sto df2.  

Zenith angle calculation formulas: 
TST=LT+ET+LC
ET=9.87sin(2B)-7.53cos(B)-1.5sin(B) , B=(360/364)(N-81)
LC=-4(lon_m-lon)
h=15(TST-12)
d=23.45sin[(360/364)(284+N)]
cozZ=sin(lat)sin(d)+cos(lat)cos(d)cos(h)
Z=arccos(cosZ) 

Solar constant is Sa=1367W/m^2 according to the WMO

@author: yiann
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytz



#%% Zenith angle calculation

# latitude (lat) and longitude (lon)
lat = 38.29136389  # (>0 north & <0 south) 
lon = 21.78858333  # (>0 east & <0 west)
# longitude of location's standard meridian (prime meridian at 0)
lon_m = 30.0  # Greece is +2 hours from 0 (15deg per standard meridian)

# inserting datetimes
data_file = 'Solar_1min_2021.txt'
df = pd.read_csv(data_file, usecols=[0], header=None, parse_dates=True, names=['Datetime'])

# df.columns = ['Datetime'] Μπορείς να το ορίσεις στην προηγούμενη εντολή

# convert 'Datetime' to datetime object (alliws den doulevei tipota) Λογικό!
# gia ton ypologismo ths zenitheias gwnias doulevw se topikh wra rologiou 
df['Datetime'] = pd.to_datetime(df['Datetime'])
# day count for the year
df['Day Count'] = df['Datetime'].dt.dayofyear


# local time (LT) calculation  
# number of hours in a day
df['Hours'] = df['Datetime'].dt.hour
# number of minutes in an hour
df['Minutes'] = df['Datetime'].dt.minute
# LT in decimal format   
df['LT'] = df['Hours'] + df['Minutes'] / 60

# equation (correction) of time in decimal format (ET Decimal)
df['beta'] = (360 / 364) * (df['Day Count'] - 81)
df['ET'] = (9.87 * np.sin(np.deg2rad(2 * df['beta'])) - 7.53 * np.cos(np.deg2rad(df['beta']))
            - 1.5 * np.sin(np.deg2rad(df['beta']))) / 60

# longitude correction (LC) calculation, in decimal hours format
df['LC'] = ((-4) * (lon_m - lon)) / 60

# true solar time (TST) calculation in decimal format
df['TST'] = df['LT'] + df['ET'] + df['LC']
# making negative values positive (alla allazei kai h stoixhsh se sxesh me thn mera?)
df['TST'] = np.where(df['TST'] < 0, df['TST'] + 24, df['TST'])

# hour angle calculation (h - in degrees)
df['h'] = (df['TST'] - 12) * 15

# declination of the sun calculation (delta - in degrees)
df['delta'] = 23.45 * np.sin(np.deg2rad((360 / 365) * (df['Day Count'] + 284)))

# cozZ calculation - where Z is the zenith angle
df['cosZ'] = np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(df['delta'])) + np.cos(np.deg2rad(lat)) * \
             np.cos(np.deg2rad(df['delta'])) * np.cos(np.deg2rad(df['h']))

# zenith angle (Z - in degrees) calculation
df['Z'] = np.rad2deg(np.arccos(df['cosZ']))

# making the index a UTC+00 Datetime (Datetime is EET/UTC+02 - local time)   
df.index = df['Datetime']
# df = df.rename_axis(index = 'Datetime UTC')
# df.index.rename('Daytime UTC', inplace=True)
df.index.name = 'Datetime UTC'
grc_std = pytz.timezone('Etc/GMT-2')  # EET (UTC+02)
utc = pytz.timezone('UCT')  # UTC+00
df.index = df.index.tz_localize(grc_std).tz_convert(utc)



#%% GHI and zenith angle values input

# irradiance values input 
df2 = pd.read_csv(data_file, index_col=[0], usecols=[0, 6], sep=',', header=None,
                  parse_dates=True, na_values='"NAN"')
df2.columns = ['GH']
df2['GH'] = df2['GH'] * 1000 / 8.63  # Converts mV to W/m^2
df2.index = df.index  # anagkaio? (des arxiko comment gia duplicate indexes)
df2['Z'] = df['Z']
df2.dropna()  # Deletes missing values

'''
to parakatw den xreiazetai an parw ta index apo to df kai ta valw sto df2:

# make the index a UTC+0 DatetimeIndex (original/local Datetime is UTC+2) 
df2.index.names = ['Datetime UTC']   
grc_std = pytz.timezone('Etc/GMT-2')     # Greek Winter Timezone (+2h offset)      
utc=pytz.timezone('UCT')     # UTC (no offset)    
df2.index = df2.index.tz_localize(grc_std).tz_convert(utc)
'''


# %% Test Plots

# aplo plotarisma twn timwn ghi   
plt.plot(df2['GH'])
plt.show()
# aplo plotarisma ths zenitheias 
plt.plot(df2['Z'])
plt.show()
# aplo plotarima ths apoklishs
plt.plot(df['delta'])
plt.show()
# aplo plotarisma ths diorthwshs xronou 
plt.plot(df['ET'])
plt.show()
