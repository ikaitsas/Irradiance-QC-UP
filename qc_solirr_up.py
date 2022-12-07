# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

This code is an attempt to perform quality control on a solar irradiance
dataset. The GH & DIFF values are imported in mV. They are converted to W/m2
in the code. The DNI values do not need conversion. The datetime column of
df (where the zenith angle is calculated) is in Greek Winter Time (UTC+02)
at first and then its index is  converted to UTC+00 datetime index.

Zenith angle calculation formulas:
TST=LT+ET+LC
ET=9.87sin(2B)-7.53cos(B)-1.5sin(B) , B=(360/364)(N-81)
LC=-4(lon_m-lon)
h=15(TST-12)
d=23.45sin[(360/364)(284+N)]
cozZ=sin(lat)sin(d)+cos(lat)cos(d)cos(h)
Z=arccos(cosZ)

The Earth-Sun distance adjusted mean solar constant is Sav=1366 W/m^2.
H timh auth einai parmenh apo tis diafaneies tou 'Fysikh Atmosfairas II'.
The formula for its variation through the year is:
Sa=Sav[1+0.033cos(360n/365)]   (apo 'Systhmata Hliakhs Energeias')

All negative GH(mV) or DIF(mV) values are set to 0. And rejected (flag==-1).
Datapoints that have Z>80 are also flagged with -1.
Datapoints with flag==-1 do not undergo any testing, as they are rejected.
Negative DN values are not flagged or corrected (this might need change).
It should be noted that DN is estimated from the sunshine duration, not
measured directly by a pyrheliometer or other device.

The data are flagged as follows:
flag            meaning
-1              Z>80 and/or original GH or DIF values negative
0               datapoint passes all QC tests
1               datapoint fails QC1
2               datapoint passes QC1 but fails QC2
3               datapoint passes QC1 & QC2 but fails QC3

Datapoints that fail QC2 are not tested any further.
It is not specified if datapoints that fail QC3 fail because of the closure
ratio test, the diffuse ratio test, or both. Might change flag values to
differentiate betwwen the 2 tests.

The test results and limit values are all stored in dataframe "df2". It is then
saved in a .txt file, as it is. The file's/dataframe's column names mean:
UTC: datetime in UTC -- corresponds to a datapoint
DIF: diffuse horizontal irradiance
GH: gobal horizontal irradiance
DN: direct normal irradiance
flag: the datapoint flag, the meaning of each value is mentioned above
Z: solar zenith angle
m0: cosine of Z (cosZ)
Day Count: Number of the day in the year -- January 1st = 1
Sa: extraterrestrial irradiance (at the top of the atmosphere). Day Count is
    used for its calculation through the year
ppl_gh: physical possible limit for GH
ppl_dif: physical possible limit for DIF
erl_gh: extremely rare limit for GH
erl_dif: extremely rare limit for DIF
sumw: DN * cosZ + DIF, calculated from the data
closr: closure ratio -- the closure equation is GH = DNcosZ+DIF, with the ratio
       here being closr=GH/(DNcosZ+DIF), calculated from the data
dif_r: diffuse ratio -- DIF/GH

In the __flagged.txt file, NaN values are stored as '', which are empty(?)
spaces. That might be changed to something like 'NaN' or 'null', I'll see. It 
needs to be specified in what string/form Nan values are stored probably, maybe 
in the _stats file?? I ideally want columns to not be objects...
Same happens in the __qc.csv file.

The complete timeseries (from XXXX-01-01 00:00:00 to XXXX-12-31 23:59:00) is
stored in a __qc.csv file. Only the timeseries, GHI & DIF values are stored in
it. Timestamps that were missing in the original data have GH, DIF inserted as
NaNs. Timestamps with flag == 1 have their GH, DIF values converted to NaNs.
Timestamps with flag == -1 have their GH, DIF values converted to 0.

@author: yiann
"""
#%% Import Stuff

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytz
from seaborn import set_style
from datetime import datetime, timedelta

set_style('whitegrid')  # vrhka to grid style tou paper lol


'''  USER INPUT - START  '''

# latitude (lat) and longitude (lon)
lat = 38.29136389  # (>0 north & <0 south)
lon = 21.78858333  # (>0 east & <0 west)

# location's hourly offset from prime merindian (h_off=0 there - negative west)
h_off = 2
local_tz = 'Etc/GMT-2'  # check pytz module for more inputs

# import the year that coresponds to the data (must be integer)
data_year = 2014

'''   USER INPUT - END   '''


# configure which columns are imported from the data file
cols_needed = [0, 4, 6, 8]  # used for read_csv: "use_cols" argument
col_datetime = [0]  # used for read_csv: "index_col" argument
col_names = ['DIF', 'GH', 'DN']  # which column is which type of irradiance

closr_test = 'YES'
if data_year < 2014:
    closr_test = 'NO'
    col_names.remove('DN')
    cols_needed.remove(8)


#%% Zenith Angle Calculation

print('\nCalculating Z & importing the data...')

# generate daterange for Z calculation
dtr = pd.date_range(start=datetime(year=data_year, month=1, day=1),
                    end=datetime(year=data_year+1, month=1, day=1) + \
                        timedelta(minutes=60*h_off-1),
                    freq='min')
# the 60*h_off-1 is written that way because the start of next year is already
# present in the end argument of pd.date_range

# create datetime dataframe for specified year
df=pd.DataFrame(dtr, columns=['Datetime'])

# zenith angle is calculated in local clock time, which is converted to UTC
# in the end, at the index, so it can be joined to the irradiance data

df['Datetime'] = pd.to_datetime(df['Datetime'])
df['Day Count'] = df['Datetime'].dt.dayofyear

# local time (LT)
df['Hours'] = df['Datetime'].dt.hour
df['Minutes'] = df['Datetime'].dt.minute
df['LT'] = df['Hours'] + df['Minutes'] / 60  # decimal format

# equation (correction) of time in decimal format (ET - decimal format)
df['beta'] = (360 / 364) * (df['Day Count'] - 81)
df['ET'] = (9.87 * np.sin(np.deg2rad(2 * df['beta'])) - 7.53 * np.cos(np.deg2rad(df['beta']))
            - 1.5 * np.sin(np.deg2rad(df['beta']))) / 60

# longitude of location's standard meridian (prime meridian at 0)
lon_m = 15.0 * h_off  # (15deg per standard meridian)
# longitude correction (LC), in decimal hours format
df['LC'] = ((-4) * (lon_m - lon)) / 60

# true solar time (TST) in decimal format
df['TST'] = df['LT'] + df['ET'] + df['LC']
# making negative values positive (for consistency's sake)
df['TST'] = np.where(df['TST'] < 0, df['TST'] + 24, df['TST'])

# hour angle (h - in degrees)
df['h'] = (df['TST'] - 12) * 15

# declination of the (delta - in degrees)
df['delta'] = 23.45 * np.sin(np.deg2rad((360 / 365) * (df['Day Count'] + 284)))

# cozZ calculation - where Z is the zenith angle
df['cosZ'] = np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(df['delta'])) + np.cos(np.deg2rad(lat)) * \
             np.cos(np.deg2rad(df['delta'])) * np.cos(np.deg2rad(df['h']))

# zenith angle (Z - in degrees)
df['Z'] = np.rad2deg(np.arccos(df['cosZ']))

# making the index a UTC+00 Datetime (Datetime is in local time - EET/UTC+02)
df.index = df['Datetime']
df.index.name = 'Datetime UTC'
df.index = df.index \
    .tz_localize(pytz.timezone(local_tz)) \
    .tz_convert(pytz.timezone('UCT'))


#%% GHI, DNI, DIF, Sa and Zenith Angle Values Input

# Irradiance types
# Global Horizontal = GH
# DIffuse Horizontal = DIF
# Direct Normal = DN

data_file = f'Solar_1min_{data_year}.txt'  # generate irradiance data file
df2 = pd.read_csv(data_file, index_col=col_datetime, usecols=cols_needed,
                  sep=',', header=None, parse_dates=True, na_values="NAN")
df2.columns = col_names

# check if df2.index is not datetimeindex, to avoid errors
if df2.index.inferred_type != 'datetime64':
    df2.index=pd.to_datetime(df2.index.astype(str), errors='coerce')
    
df2.index = df2.index.tz_localize('UTC')
df2['flag'] = -1  # initialising test flags - reject all at first (flag==-1)
df2.loc[(df2['GH']>0) & (df2['DIF']>0), 'flag'] = 0  # reject GH<0 or DIF<0
df2.loc[df2['GH']<0, 'GH'] = 0  # make negative voltage values zero
df2.loc[df2['DIF']<0, 'DIF'] = 0
# maybe this has to be done to DN too??
df2['DIF'] = df2['DIF'] * 1000 / 8.64  # Converts mV to W/m^2
df2['GH'] = df2['GH'] * 1000 / 8.63

df2 = df2.join(df['Z'])  # left join is default
df2['m0'] = np.cos(np.deg2rad(df2['Z'])).clip(lower=0)  # m0=cosZ>=0
df2.loc[df2['Z']>80, 'flag'] = -1  # flag == -1 for Z>80

# drop duplicate datetimes (pref. after join operations)
df2=df2[~df2.index.duplicated(keep='first')]

# solar constant calculation (Sav = 1366 W/m2)
df2['Day Count'] = df2.index.dayofyear
df2['Sa'] = 1366 * ( 1 + 0.033 * np.cos( np.deg2rad((360 * df2['Day Count'])/365)) )
# or it can be set constant: (Sa=1366 W/m^2): df2['Sa'] = 1366

df2 = df2[df2.index.notnull()]  # deletes NaT indices
df2 = df2.dropna(subset=col_names)  # Deletes NaN/missing values


#%% Limit Tests (QC1 & QC2)

print('Preparing QC1, QC2 & QC3 figures...')

# PPL tests (QC1)
df2['ppl_gh'] = df2['Sa'] * 1.5 * df2['m0']**1.2 + 100
df2['ppl_dif'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 50
# ERL tests (QC2)
df2['erl_gh'] = df2['Sa'] * 1.2 * df2['m0']**1.2 + 50
df2['erl_dif'] = df2['Sa'] * 0.75 * df2['m0']**1.2 + 30

# plotting QC1 & QC2 (only values that correspond to Z<80 are shown)
# GH plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_gh'][df2['Z']<80], s=0.002, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_gh'][df2['Z']<80], s=0.002)
plt.scatter(df2['Z'][df2['Z']<80], df2['GH'][df2['Z']<80], s=0.002, c='k')
plt.ylim([0,2125])
plt.title('Global Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=150, loc='upper right')
plt.savefig('qc1_2_ghi_vs_zenith__'+data_file[:-4]+'.png')
plt.show()

# DIF plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_dif'][df2['Z']<80], s=0.002, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.002)
plt.scatter(df2['Z'][df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.002, c='k')
plt.ylim([0,1300])
plt.title('Diffuse Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=150, loc='upper right')
plt.savefig('qc1_2_dif_vs_zenith__'+data_file[:-4]+'.png')
plt.show()

# flag data that fail PPL test with 1
df2.loc[(df2['flag']==0) & (df2['GH']>df2['ppl_gh']), 'flag'] = 1
df2.loc[(df2['flag']==0) & (df2['DIF']>df2['ppl_dif']), 'flag'] = 1

# flag data that pass PPL test but fail ERL test with 2
df2.loc[(df2['flag']==0) & (df2['GH']>df2['erl_gh']), 'flag'] = 2
df2.loc[(df2['flag']==0) & (df2['DIF']>df2['erl_dif']), 'flag'] = 2


#%% Comparison Tests (QC3)

# closure equation test
if closr_test == 'YES':
    df2['sumw'] = df2['DN'] * df2['m0'] + df2['DIF']
    df2['closr'] = df2['GH'] / df2['sumw']
    df2.loc[(df2['Z']>80) | (df2['sumw']<50), 'closr'] = np.nan
    # sto ena paper exei ghi>50 kai sto allo sumw>50, ti na valw? tha rwthsw..

    plt.scatter(df2['Z'], df2['closr'], s=0.002, c='k')
    # ratio limits visualization
    plt.hlines(y=1.08, xmin=10, xmax=75, color='r')
    plt.hlines(y=1.15, xmin=75, xmax=85, color='r')
    plt.vlines(x=75, ymin=1.08, ymax=1.15, color='r')
    plt.hlines(y=0.92, xmin=10, xmax=75, color='r')
    plt.hlines(y=0.85, xmin=75, xmax=85, color='r')
    plt.vlines(x=75, ymin=0.85, ymax=0.92, color='r')
    plt.ylim([0,2])  # so that extreme outliers dont affect the figure..
    plt.title('Closure Ratio (GH/SUMW)')
    plt.xlabel('Zenith Angle [°]')
    plt.legend(['Data', 'Limits'], markerscale=150, loc='upper right')
    plt.savefig('qc3_ghi_sumw_ratio__'+data_file[:-4]+'.png')
    plt.show()
    
    # flag data that fail closure ratio test with 3
    df2.loc[(df2['flag']==0) & (df2['Z']<75) & (df2['closr']<0.92), 'flag'] = 3
    df2.loc[(df2['flag']==0) & (df2['Z']<75) & (df2['closr']>1.08), 'flag'] = 3
    df2.loc[(df2['flag']==0) & (df2['Z']>75) & (df2['closr']<0.85), 'flag'] = 3
    df2.loc[(df2['flag']==0) & (df2['Z']>75) & (df2['closr']>1.15), 'flag'] = 3
else:
    df2['closr'] = np.nan


'''
# flag data that fail closure ratio test with 3
df2.loc[(df2['flag']==0) & (df2['Z']<75) &
        ((df2['closr']<0.92) | (df2['closr']>1.08)), 'flag'] = 3
df2.loc[(df2['flag']==0) & (df2['Z']>75) &
        ((df2['closr']<0.85) | (df2['closr']>1.15)), 'flag'] = 3
'''
# diffuse ratio test
df2['dif_r'] = df2['DIF'] / df2['GH']
df2.loc[(df2['Z']>80) | (df2['GH']<50), 'dif_r'] = np.nan

plt.scatter(df2['Z'], df2['dif_r'], s=0.002, c='k')
# ratio limits visualization
plt.hlines(y=1.05, xmin=10, xmax=75, color='r')
plt.hlines(y=1.1, xmin=75, xmax=85, color='r')
plt.vlines(x=75, ymin=1.05, ymax=1.1, color='r')
plt.ylim([0,1.4])
plt.title('Diffuse Ratio (DIF/GH)')
plt.xlabel('Zenith Angle [°]')
plt.savefig('qc3_dif_ghi_ratio__'+data_file[:-4]+'.png')
plt.show()

# flag remaining data that fail the diffuse ratio test with 3
df2.loc[(df2['flag']==0) & (df2['Z']<75) & (df2['dif_r']>1.05), 'flag'] = 3
df2.loc[(df2['flag']==0) & (df2['Z']>75) & (df2['dif_r']>1.10), 'flag'] = 3


# %% Timeseries Plots

print('\nPlotting the timeseries of the data...')

# plotting of irradiance data and limits versus time/date
plt.scatter(df2.index[df2['Z']<80], df2['ppl_gh'][df2['Z']<80],
            s=0.001, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_gh'][df2['Z']<80], s=0.001)
plt.scatter(df2.index[df2['Z']<80], df2['GH'][df2['Z']<80], s=0.005, c='k')
plt.ylim([0,2125])
plt.title('GHI')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=80, loc='upper right')
plt.savefig('ghi_vs_time__'+data_file[:-4]+'.png')
plt.show()

#plt.hist(df2['GH'][df2['Z']<80], bins=60, log=True, color='k')
#plt.show()

if closr_test == 'YES':
    plt.scatter(df2.index[df2['Z']<80], df2['DN'][df2['Z']<80], s=0.02, c='k')
    plt.title('DNI')
    plt.ylabel('Irradiance [$W/m^2$]')
    #plt.savefig('dni_vs_time__'+data_file[:-4]+'.png')
    plt.show()
    
    #plt.hist(df2['DN'][df2['Z']<80], bins=50, log=True, color='k')
    #plt.show()

plt.scatter(df2.index[df2['Z']<80], df2['ppl_dif'][df2['Z']<80],
            s=0.001, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.001)
plt.scatter(df2.index[df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.005, c='k')
plt.ylim([0,1300])
plt.title('DIF')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=80, loc='upper right')
plt.savefig('dif_vs_time__'+data_file[:-4]+'.png')
plt.show()

#plt.hist(df2['DIF'][df2['Z']<80], bins=50, log=True, color='k')
#plt.show()


#%% Results Exportation

# display flagging results in the console
print('\nSome pass/fail stats for', data_file[:-4] + ':',
      '\n\nTotal number of datapoints:',
      str(len(df2.index)),
      '\nTotal number of non-eligible datapoints:',
      str(len(df2[df2['flag']==-1].index)),
      '\nTotal number of tested datapoints:',
      str(len(df2[df2['flag']>-1].index)),
      '\n\nDatapoints that fail test 1:',
      str(len(df2[df2['flag']==1].index)),
      '\n\nDatapoints that pass test 1 but fail test 2:',
      str(len(df2[df2['flag']==2].index)),
      '\n\nDatapoints that pass tests 1 & 2 but fail test 3:',
      str(len(df2[df2['flag']==3].index)),
      '\nOf those that fail test 3, these many fail in:',
      '\n\tThe closure ratio test:',
      str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                  ((df2['closr']<0.92) | (df2['closr']>1.08))].index)+ \
          len(df2[(df2['flag']==3) & (df2['Z']>75) &
                  ((df2['closr']<0.85) | (df2['closr']>1.15))].index)),
      '\n\tThe diffuse ratio test:',
      str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                  (df2['dif_r']>1.05)].index)+ \
          len(df2[(df2['flag']==3) & (df2['Z']>75) &
                  (df2['dif_r']>1.1)].index))
      )

print('\n\nPlease wait for the exportation of the data & results.',
      '\nA "Done" message will pop up in the terminal upon completion.',
      '\nIt might take some time...'
      )

# export the pass\fail stats about the data to a text file
with open(data_file[:-4]+'__flagged_stats.txt', 'w') as yearly_stats_file:
    print('\n\nSome pass/fail stats for', data_file[:-4] + ':',
          '\n\nTotal number of datapoints:',
          str(len(df2.index)),
          '\nTotal number of non-eligible datapoints:',
          str(len(df2[df2['flag']==-1].index)),
          '\nTotal number of tested datapoints:',
          str(len(df2[df2['flag']>-1].index)),
          '\n\nDatapoints that fail test 1:',
          str(len(df2[df2['flag']==1].index)),
          '\n\nDatapoints that pass test 1 but fail test 2:',
          str(len(df2[df2['flag']==2].index)),
          '\n\nDatapoints that pass tests 1 & 2 but fail test 3:',
          str(len(df2[df2['flag']==3].index)),
          '\nOf those that fail test 3, these many fail in:',
          '\n\tThe closure ratio test:',
          str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                      ((df2['closr']<0.92) | (df2['closr']>1.08))].index)+ \
              len(df2[(df2['flag']==3) & (df2['Z']>75) &
                      ((df2['closr']<0.85) | (df2['closr']>1.15))].index)),
          '\n\tThe diffuse ratio test:',
          str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                      (df2['dif_r']>1.05)].index)+ \
              len(df2[(df2['flag']==3) & (df2['Z']>75) &
                      (df2['dif_r']>1.1)].index)),
          file=yearly_stats_file
          )
        
with open(data_file[:-8]+'__all_years_summary.txt', 'a') as stats_file:
    print('\n\nSome pass/fail stats for', data_file[:-4] + ':',
          '\n\nTotal number of datapoints:',
          str(len(df2.index)),
          '\nTotal number of non-eligible datapoints:',
          str(len(df2[df2['flag']==-1].index)),
          '\nTotal number of tested datapoints:',
          str(len(df2[df2['flag']>-1].index)),
          '\n\nDatapoints that fail test 1:',
          str(len(df2[df2['flag']==1].index)),
          '\n\nDatapoints that pass test 1 but fail test 2:',
          str(len(df2[df2['flag']==2].index)),
          '\n\nDatapoints that pass tests 1 & 2 but fail test 3:',
          str(len(df2[df2['flag']==3].index)),
          '\nOf those that fail test 3, these many fail in:',
          '\n\tThe closure ratio test:',
          str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                      ((df2['closr']<0.92) | (df2['closr']>1.08))].index)+ \
              len(df2[(df2['flag']==3) & (df2['Z']>75) &
                      ((df2['closr']<0.85) | (df2['closr']>1.15))].index)),
          '\n\tThe diffuse ratio test:',
          str(len(df2[(df2['flag']==3) & (df2['Z']<75) &
                      (df2['dif_r']>1.05)].index)+ \
              len(df2[(df2['flag']==3) & (df2['Z']>75) &
                      (df2['dif_r']>1.1)].index)),
          file=stats_file
          )

# export the data to a text file
df2.index=df2.index.tz_localize(None)
df2=df2.round(decimals=3)
#pd.set_option('display.precision', 3)
df2.to_csv(data_file[:-4]+'__flagged.txt', index_label='UTC')

# export a chronologicaly complete timeseries of the data
# needs polishing for if there is no XXXX-12-31 23:59:00
# all XXXX-01-01 00:00:00 to XXXX-21-31 23:59:00 timestamps must be contained
# maybe with joining with a datataframe/datetimeindex?
# but that will make the code pretty 'heavy' (it already is..)

if df2.index[0] != (pd.to_datetime(f'{data_year}-1-1 00:00:00') or \
                    pd.to_datetime(f'{data_year}-12-31 23:59:00')):
    df3 = df2[['GH','DIF','flag']] \
    .reindex(pd.date_range(f'{data_year}-1-1 00:00:00',
                           f'{data_year}-12-31 23:59:00',
                           freq='1min'), fill_value=np.nan)
else:
    df3 = df2[['GH','DIF','flag']].asfreq(freq='1min')

df3.loc[df3['flag'] == -1, 'GH'] = 0
df3.loc[df3['flag'] == -1, 'DIF'] = 0
df3.loc[df3['flag'] == 1, 'GH'] = np.nan
df3.loc[df3['flag'] == 1, 'DIF'] = np.nan
df3 = df3.drop(columns=['flag'])
df3.to_csv(data_file[:-4]+'__qc.csv', index_label='UTC')

with open(data_file[:-4]+'__flagged_stats.txt', 'a') as yearly_stats_file:
    print('\nNumber of timestamps with missing values:',
          str(len(df3.index)-len(df2.index)),
          file=yearly_stats_file
          )
    
with open(data_file[:-8]+'__all_years_summary.txt', 'a') as stats_file:
    print('\nNumber of timestamps with missing values:',
          str(len(df3.index)-len(df2.index)),
          file=stats_file
          )

print('\n\nDone.')    
'''
WARNING
rounding changes a bit the flagging & test results. Need to fix this.

issue: rounding seems to change a bit the number of datapoints that pass/fail
each test
e.g. without rounding, 100990 fail closure ratio test
     with 3 decimals rounding, 100484
     with 4 decimals rounding, 100947
     with 5 decimals rounding, 100989
     
the tests and stats are all done/taken with unrounded data in this code, but
they are stored in the .csv file rounded
running the tests/stats with rounded data produces slightly different stats

this might be happening because i have put only the "<" or ">" condition, not
the "=<" or ">=" condition for the tests, so the rounding towards the limits
makes a datapoint evade flagging.
'''