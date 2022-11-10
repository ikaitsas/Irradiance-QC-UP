# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

O kwdikas kanei elegxo twn ghi timwn
H ghi antistoixei sthn sthlh 6 apo tis 0-13 (14 sthles to arxeio)
Den emfanizetai arxika h ghi, alla h tash tou pyranometrou

ΣΥΝΕΧΙΖΟΥΜΕ ΣΤΟΝ ΓΗΙ_ΤΕΣΤΣ2

@author: yiann
"""

import numpy as np
import pandas  as pd
import matplotlib.pyplot as plt 

# eisagwgh tou arxeiou ston kwdika
data_file='Solar_1min_2021copy.txt'
old_data=pd.read_csv(data_file,index_col=[0],usecols=[0,6],header=None,parse_dates=[0])
old_data.columns=['mv'] #


# symplhrwsh twn leptwn pou leipoun me nan
filled_data=old_data.asfreq('1min')

# metatroph twn mv se monades entashs
filled_data['ghi']=(1000/8.63)*filled_data['mv']


# aplo plotarisma twn timwn ghi
filled_data.plot(y='ghi') 














