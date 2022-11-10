# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

O kwdikas kanei elegxo twn ghi timwn
H ghi antistoixei sthn sthlh 6 apo tis 0-13 (14 sthles to arxeio)
Den emfanizetai arxika h ghi, alla h tash tou pyranometrou

@author: yiann
"""
import pandas  as pd
import matplotlib.pyplot as plt 

# eisagwgh tou arxeiou ston kwdika
data_file='Solar_1min_2021.txt'
df=pd.read_csv(data_file, index_col=[0], usecols=[0,4,6,8], sep=',', header=None, 
    parse_dates=True, na_values='"NAN"')
df.columns=['DIF','GH','DN']
df.dropna() # Deletes missing values
#df = df * 1000 / 8.63 # Converts mV to W/m^2

# symplhrwsh twn leptwn pou leipoun me nan
# filled_data=old_data.asfreq('1min') Do not touch missing data


# aplo plotarisma twn timwn ghi  
plt.plot(df['GH'])
plt.show()

plt.plot(df['DN'])
plt.show()

plt.plot(df['DIF'])
plt.show()












