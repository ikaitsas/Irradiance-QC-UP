# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

This code is an attempt to perform quality control in a solar irradiance 
dataset. The GH & DIFF values are imported in mV. They are converted to W/m2
in the code. The DN values do not need conversion. The datetime column/index 
of df (where the zenith angle is calculated) is in Greek Winter Time (UTC+2). 
 

Zenith angle calculation formulas: 
TST=LT+ET+LC
ET=9.87sin(2B)-7.53cos(B)-1.5sin(B) , B=(360/364)(N-81)
LC=-4(lon_m-lon)
h=15(TST-12)
d=23.45sin[(360/364)(284+N)]
cozZ=sin(lat)sin(d)+cos(lat)cos(d)cos(h)
Z=arccos(cosZ) 


** added 2022/10/19 ** 

The Earth-Sun distance adjusted mean solar constant is Sav=1366 W/m^2.
H timh auth einai parmenh apo tis diafaneies tou 'Fysikh Atmosfairas II'. 
The formula for its variation through the year is: 
Sa=Sav[1+0.033cos(360n/365)] , n=1 for Jan 1st. 

Ta diagrammata tha ta peripoihthw sto telos. 

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

# inserting datetimes from the irradiance file 
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

# thn metatroph apo UTC+2 se UTC na thn kanw afou teleiwsw tous elegxous?? 


'''  apo edw kai katw yparxoun prosthikes ''' 
#%% GHI, DNI, DIF, Sa and zenith angle values input


# irradiance values input  
df2 = pd.read_csv(data_file, index_col=[0], usecols=[0, 4, 6, 8], sep=',', header=None,
                  parse_dates=True, na_values='"NAN"')
df2.columns = ['DIF', 'GH', 'DN']  # diffuse-global horizontal-direct normal
df2['DIF'] = df2['DIF'] * 1000 / 8.64  # Converts mV to W/m^2 
df2['GH'] = df2['GH'] * 1000 / 8.63  # Converts mV to W/m^2
df2.index = df.index  # anagkaio (des arxiko comment gia duplicate indexes)
df2['Z'] = df['Z']  # isws xreiazetai na eisagw kai to cosZ edw?? 
# it must be m0(=cosZ)>0 (or Z<90) for the QC tests 
df2['m0'] = df['cosZ'].clip(lower=0) 

# solar constant calculation (Sav = 1366 W/m2)
df2['Sa'] = 1366 * ( 1 + 0.033 * np.cos( np.deg2rad((360 * df['Day Count'])/365)) ) 
# alliws an thn thewrhsw statherh Sa=1366 W/m^2:
#df2['Sa'] = 1366 
df2.dropna()  # Deletes missing values 

'''
Κοιταξα μερικες τιμες και αρκετες ειναι πολυ εκτος σε σχεση με το PPL, κωριως 
στις ωρες που ανατελλει ο ηλιος. Επισης, οι τιμες πεφτουν πολυ νωριτερα σε
σχεση με αυτο. Η ασθητη αποκλιση ειναι περι τις 2 ωρες. Μαλλον δεν εχω κανει
καλη αντιστοιχια με το δφ, οπου υπολογισα την ζενιθεια γωνια...

Στο αρχικο αρχειο η ωρα ειναι σε ΘΤΨ. Οποτε και τ οδφ2 ειναι σε ΘΤΨ. Ο αληθης
ηλιακος χρονος υπολογιζεται με βαση την διορθωση της τοπικης ωρας απο τους 
τυπους που χρησιμοποιω. Ακομα και αν εισαγω τις αρχικες ημερομηνιες απο το 
αρχειο, οι οποιες ειναι σε ΘΤΨ, οι υπολογισμοι τη αναγνωριζουν ως τοπικη,
δηλαδη ως ΘΤΨ+2. Παρολο που εισηγαγα τις ωρες απο την ιδια πηγη, το δφ την
αντιλαμβανεται ως τοπικη, που δεν ειναι φυσικα στο αρχικο αρχειο, αρα και στο
δφ2. Οποτε θελει μετατροπη στο τελος. Το εχω κανει.

Η βλακεια που νομιζω εχω κανει ειναι η επεγγραφη/αντικατασταση της αρχικης
ωρας, ως ΘΤΨ στο δφ2, με αυτην του δφ, οπου μετα την μετατροπη αρχιζει απο την
2020/12/31 22:00:00, 2 ωρες νωριτερα δηλαδη. Αρα και η ζενιθεια γωνια, που 
εισαγεται απο το δφ, παει ξανα 2 ωρες πισω. Για αυτο παρατηρειται η διαφορα 
των 2 ωρων, που δεν θα επρεπε να υπαρχει. 

Πρεπει να βρω εναν τροπο να αντιστοιχησω σωστα στην ζενιθεια γωνια με τις 
τιμες της εντασης της ακτινοβολιας. 
''' 



#%% Limits test 


# PPL tests (QC1) 
df2['ppl_min'] = -4  # minimun value for PPL test - common to all 
df2['ppl_gh'] = df2['Sa'] * 1.5 * df2['m0']**1.2 + 100  # max PPL for GH
df2['ppl_dif'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 50  # max PPL for DIF
df2['ppl_dn'] = df2['Sa']  # max PPL for DN 

# ERL tests (QC2) 
df2['erl_min'] = -2  # minimun value for ERL test - common to all 
df2['erl_gh'] = df2['Sa'] * 1.2 * df2['m0']**1.2 + 50  # max ERL for GH 
df2['erl_dif'] = df2['Sa'] * 0.75 * df2['m0']**1.2 + 30  # max ERL for DIF 
df2['erl_dn'] = df2['Sa'] * 0.95 * df2['m0']**0.2 + 10  # max ERL for DN 
 
# mporw na kanw kai gia to DNcosZ (amesh katheth aktinobolia) an xreiazetai 

# plotting QC1 & QC2 

# GH plots 
# only values that correspond to Z<90 (m0>0) are shown 
plt.scatter(df2['Z'][df2['Z']<90], df2['GH'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['ppl_gh'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['erl_gh'][df2['Z']<90], s=1) 
plt.show() 

# DIF plots 
# only values that correspond to Z<90 (m0>0) are shown 
plt.scatter(df2['Z'][df2['Z']<90], df2['DIF'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['ppl_dif'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['erl_dif'][df2['Z']<90], s=1) 
plt.show() 

# DN plots 
# only values that correspond to Z<90 (m0>0) are shown 
plt.scatter(df2['Z'][df2['Z']<90], df2['DN'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['ppl_dn'][df2['Z']<90], s=1) 
plt.scatter(df2['Z'][df2['Z']<90], df2['erl_dn'][df2['Z']<90], s=1) 
plt.show() 



#%% Comparison tests 

#




# %% Test Plots


# aplo plotarisma ths zenitheias 
plt.plot(df2['Z'])
plt.show()
# aplo plotarima ths apoklishs
plt.plot(df['delta'])
plt.show()
# aplo plotarisma ths diorthwshs xronou 
plt.plot(df['ET'])
plt.show()
# aplo plotarisma twn timwn aktinovolias    
plt.plot(df2['GH'])
plt.show()
plt.plot(df2['DN'])
plt.show()
plt.plot(df2['DIF'])  
plt.show() 
