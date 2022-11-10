# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 17:53:18 2022

@author: yiann
"""

import pandas  as pd 
import pytz

# eisagwgh tou arxeiou ston kwdika
data_file='Solar_1min_2021.txt'
df=pd.read_csv(data_file, index_col=[0], usecols=[0,6], sep=',', header=None, 
    parse_dates=True, na_values='"NAN"')
df.columns=['GH']
df.dropna() # Deletes missing values
df = df * 1000 / 8.63 # Converts mV to W/m^2

# dokimazw timezones apo to pytz
grc = pytz.timezone('Europe/Athens') # includes DST, den mou kanei
grc_w = pytz.timezone('Etc/GMT-2') # auto pou psaxnw???
utc=pytz.timezone('UCT')

df.index = df.index.tz_localize(grc_w).tz_convert(utc)  

''' ------------------------ TA KATAFERAAAAAAAA? ---------------------------'''


