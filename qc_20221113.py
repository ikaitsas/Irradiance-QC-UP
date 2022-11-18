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
Datapoints thathave Z>80 are also flagged with -1.
Datapoints with flag==-1 do not undergo any testing.
Negative DN values are not flagged or corrected.
It should be noted that DN is estimated from the sunshine duration, not
measured directly by a pyrheliometer or other device.

The data are flagged as follows:
flag            meaning
-1              Z>80
0               datapoint passes all QC tests
1               datapoint fails QC1
2               datapoint passes QC1 but fails QC2
3               datapoint passes QC1 & QC2 but fails QC3

Datapoints that fail QC2 are not tested any further.
It is not specified if datapoints that fail QC3 fail because of the closure
ratio test, the diffuse ratio test, or both.

@author: yiann
"""
#%% Import Stuff


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytz
import seaborn as sns

sns.set_style('whitegrid')  # vrhka to grid style tou paper lol

# latitude (lat) and longitude (lon)
lat = 38.29136389  # (>0 north & <0 south)
lon = 21.78858333  # (>0 east & <0 west)
# longitude of location's standard meridian (prime meridian at 0)
lon_m = 30.0  # Greece is +2 hours from 0 (15deg per standard meridian)

# import irradiance data and create datetimes for zenith angle calculation
data_file = 'Solar_1min_2021.txt'
dtr=pd.date_range(start='2021-1-1 00:00:00', end='2022-1-1 1:59:00',
                  freq='min')


#%% Zenith Angle Calculation

print('\nCalculating Z & importing the data...')

# create datetime dataframe for specified year
df=pd.DataFrame(dtr, columns=['Datetime'])

# zenith angle is calculated in local clock time, which is converted to UTC
# in the end, at the index, so it can be joined to the irradiance data

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

# making the index a UTC+00 Datetime (Datetime is in local time - EET/UTC+02)
df.index = df['Datetime']
df.index.name = 'Datetime UTC'
grc_std = pytz.timezone('Etc/GMT-2')  # EET (UTC+02)
utc = pytz.timezone('UCT')  # UTC+00
df.index = df.index.tz_localize(grc_std).tz_convert(utc)


#%% GHI, DNI, DIF, Sa and Zenith Angle Values Input

# Global Horizontal = GH
# DIffuse Horizontal = DIF
# Direct Normal = DN
# Beam Irradiance = BI (BI=DNcosZ)

# irradiance values input
df2 = pd.read_csv(data_file, index_col=[0], usecols=[0, 4, 6, 8], sep=',',
                  header=None, parse_dates=True, na_values='"NAN"')
df2.columns = ['DIF', 'GH', 'DN']
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

df2 = df2.dropna(subset=['DIF', 'GH', 'DN'])  # Deletes NaN/missing values


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
plt.title('Global Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=150, loc='upper right')
plt.savefig('qc1_2_ghi_vs_zenith.png')
plt.show()

# DIF plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_dif'][df2['Z']<80], s=0.002, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.002)
plt.scatter(df2['Z'][df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.002, c='k')
plt.title('Diffuse Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=150, loc='upper right')
plt.savefig('qc1_2_dif_vs_zenith.png')
plt.show()

# flag data that fail PPL test with 1
df2.loc[(df2['flag']==0) & (df2['GH']>df2['ppl_gh']), 'flag'] = 1
df2.loc[(df2['flag']==0) & (df2['DIF']>df2['ppl_dif']), 'flag'] = 1

# flag data that pass PPL test but fail ERL test with 2
df2.loc[(df2['flag']==0) & (df2['GH']>df2['erl_gh']), 'flag'] = 2
df2.loc[(df2['flag']==0) & (df2['DIF']>df2['erl_dif']), 'flag'] = 2


#%% Comparison Tests (QC3)

# closure equation test
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
plt.ylim([0,2])  # so that extreme outliers dont affect the diagram..
plt.title('Closure Ratio (GH/SUMW)')
plt.xlabel('Zenith Angle [°]')
plt.legend(['Data', 'Limits'], markerscale=150, loc='upper right')
plt.savefig('qc3_ghi_sumw_ratio.png')
plt.show()

# flag data that fail closure ratio test with 3
df2.loc[(df2['flag']==0) & (df2['Z']<75) &
        ((df2['closr']<0.92) | (df2['closr']>1.08)), 'flag'] = 3
df2.loc[(df2['flag']==0) & (df2['Z']>75) &
        ((df2['closr']<0.85) | (df2['closr']>1.15)), 'flag'] = 3

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
plt.savefig('qc3_dif_ghi_ratio.png')
plt.show()

# flag remaining data that fail the diffuse ratio test with 3
df2.loc[(df2['flag']==0) & (df2['Z']<75) & (df2['dif_r']>1.05), 'flag'] = 3
df2.loc[(df2['flag']==0) & (df2['Z']>75) & (df2['dif_r']>1.1), 'flag'] = 3


#%% Climatological Limits Tests (QC4)

print('Preparing QC4 figures...')

# GH rejection rate visualization
reject_gh=[]
i_list=[]
for i in range(120, -5, -5):
    df2['cll_gh'] = df2['Sa'] * (i/100) * df2['m0']**1.2 + 50
    reject_gh.append(100*len(df2['GH'][(df2['GH']>df2['cll_gh'])&(df2['flag']!=-1)])/ \
                     len(df2.index[df2['flag']!=-1]))
    i_list.append(i/100)
    
plt.scatter(i_list, reject_gh, marker='+', c='k')
#plt.yscale('log')
plt.title('Rejection Rate for GHI')
plt.ylabel('Rejection Percentage')
plt.xlabel('Climatological Limits Coefficient - GHI')
plt.show()

reject_gh = np.transpose(np.array([i_list, reject_gh]))

# DIF rejection rate visualization
reject_dif=[]
j_list=[]
for j in range(75, -5, -5):
    df2['cll_dif'] = df2['Sa'] * (j/100) * df2['m0']**1.2 + 30
    reject_dif.append(100*len(df2['DIF'][(df2['DIF']>df2['cll_dif'])&(df2['flag']!=-1)])/ \
                      len(df2.index[df2['flag']!=-1]))
    j_list.append(j/100)
    
plt.scatter(j_list, reject_dif, marker='+', c='k')
#plt.yscale('log')
plt.title('Rejection Rate for DIF')
plt.ylabel('Rejection Percentage')
plt.xlabel('Climatological Limits Coefficient - DIF')
plt.show()

reject_dif = np.transpose(np.array([j_list, reject_dif]))

# ta ekana epi 100 ta coefficients giati den evgaza akrh me float values sto
# for loop, to idio pragma einai praktika kiolas
# mporei na souloupwthei, kanontas tis kenes listes kena arrays eksarxhs??


# %% Test Plots

print('Plotting the timeseries of the data...')

# plotting of irradiance data and limits versus time/date
plt.scatter(df2.index[df2['Z']<80], df2['ppl_gh'][df2['Z']<80],
            s=0.001, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_gh'][df2['Z']<80], s=0.001)
plt.scatter(df2.index[df2['Z']<80], df2['GH'][df2['Z']<80], s=0.005, c='k')
plt.title('GHI')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=80, loc='upper right')
plt.savefig('ghi_vs_time.png')
plt.show()

#plt.hist(df2['GH'][df2['Z']<80], bins=60, log=True, color='k')
#plt.show()

plt.scatter(df2.index[df2['Z']<80], df2['DN'][df2['Z']<80], s=0.02, c='k')
plt.title('DNI')
plt.ylabel('Irradiance [$W/m^2$]')
plt.savefig('dni_vs_time.png')
plt.show()

#plt.hist(df2['DN'][df2['Z']<80], bins=50, log=True, color='k')
#plt.show()

plt.scatter(df2.index[df2['Z']<80], df2['ppl_dif'][df2['Z']<80],
            s=0.001, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.001)
plt.scatter(df2.index[df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.005, c='k')
plt.title('DIF')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Data'], markerscale=80, loc='upper right')
plt.savefig('dif_vs_time.png')
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
      '\nDatapoints that pass test 1 but fail test 2:',
      str(len(df2[df2['flag']==2].index)),
      '\nDatapoints that pass tests 1 & 2 but fail test 3:',
      str(len(df2[df2['flag']==3].index))
      )

print('\n\nPlease wait for the exportation of the data & results.',
      '\nA "Done" message will pop up in the terminal upon completion.',
      '\nIt might take some time...'
      )

# export the data to a text file
df2.to_csv(data_file[:-4]+'_flagged.txt', index_label='UTC')

# export the pass\fail stats about the data to a text file
with open(data_file[:-4]+'_flagged_stats.txt', 'w') as stats_text:
    print('\nSome pass/fail stats for', data_file[:-4] + ':',
          '\n\nTotal number of datapoints:',
          str(len(df2.index)),
          '\nTotal number of non-eligible datapoints:',
          str(len(df2[df2['flag']==-1].index)),
          '\nTotal number of tested datapoints:',
          str(len(df2[df2['flag']>-1].index)),
          '\n\nDatapoints that fail test 1:',
          str(len(df2[df2['flag']==1].index)),
          '\nDatapoints that pass test 1 but fail test 2:',
          str(len(df2[df2['flag']==2].index)),
          '\nDatapoints that pass tests 1 & 2 but fail test 3:',
          str(len(df2[df2['flag']==3].index)),
          file=stats_text
          )

print('\n\nDone.')


#%% ti menei - sxolia

'''
# aplo plotarisma ths zenitheias
plt.plot(df2['Z'])
plt.show()
# aplo plotarima ths apoklishs
plt.plot(df['delta'])
plt.show()
# aplo plotarisma ths diorthwshs xronou
plt.plot(df['ET'])
plt.show()

poly o,ti na nai vgainei to ghi/sumw diagramma..
prepei na kanw kai ta pass/fail pososta na exw plhrh eikona..
mhpws epeidh einai pio 'mikra' ta DN/BI sta mikra Z apo to 'anamenomeno'..
san na fainetai na kanoun mia voutia sta mikra Z??
den paizei na peiraksa kati elpizw..
p.x. df2['DN'][df2['Z']<kapoia timh].median() goes up for Z going up
mallon vlakeies grafw..
to oti einai estimate paizei kapoio rolo??

** added 2022/10/20 **

ta duplicated indexes exoun proelthei kata to join tou Z sto df2
ta afairw afou kanw ta join
me ena aplo: df2[~df2.index.duplicated()] sto df2 eimai logika kalymmenos?

to: df2[df2.index.duplicated()] dinei ta duplicate rows
ta duplicated values kai to count tous vrhka pws dinontai apo:
count_series=df2.pivot_table(columns=['Column Name'], aggfunc='size')
count_series[count_series>1]
me index ths count_series thn sthlh me onoma 'Column Name' (ta values)
tha yparxei kai pio apodotikos tropos, alla auton skefthka twra..

** added 2022-10-22 **

df2.index.size=525561 (gia 2021), enw exw 525600 lepta se mh disekto xrono
apo to: df2['Day Count'].value_counts() vlepw apo pou leipoun metrhseis, afou
exw 24x60=1440 lepta mesa se mia hmera
p.x. leipoun 39 metrhseis sto 2021 arxeio

to: df2.drop(df2.loc[df2.index > '2021-12-31 21:59:00'].index, inplace=True)
diwxnei ta teleutaia datetimes pou einai NaN logw eisagwghs ths Z sto df2
to: df2.dropna() diwxnei ola ta pithana NaN stoixeia se kathe sthlh
df.isnull().values.any() na tsekarw an egine patata edw..

vriskw ta datetime indices opou exw NaN times ws exhs:
index = df2['Z'].index[df2['Z'].apply(np.isnan)]

** added 2022-10-29 **
mhpws na kanw ton ypologismo ths zenitheias gwnias synarthsh? tha einai kai
pio voliko gia thn enswmatwsh twn timwn ths sto df2..

h DN einai ektimwmenh sta dedomena, den ypokeitai stous alegxous autous

df2['BI'] = df2['DN'] * df2['m0']  # direct horizontal
df2['ppl_dn'] = df2['Sa']
df2['ppl_bi'] = df2['Sa'] * df2['m0']
df2['erl_dn'] = df2['Sa'] * 0.95 * df2['m0']**0.2 + 10
df2['erl_bi'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 10
# DN plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_dn'][df2['Z']<80], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_dn'][df2['Z']<80], s=0.001)
plt.scatter(df2['Z'][df2['Z']<80], df2['DN'][df2['Z']<80], s=0.001, c='k')
plt.title('DNI')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.show()
# BI plots  # mporei na mhn xreiazontai??
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_bi'][df2['Z']<80], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_bi'][df2['Z']<80], s=0.001)
plt.scatter(df2['Z'][df2['Z']<80], df2['BI'][df2['Z']<80], s=0.001, c='k')
plt.title('BI (=DNcosZ)')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.show()

** added 2022-10-30 **

df2.loc[((df2['flag']>-1) & (df2['flag']<1)) & (df2['GH']>df2['erl_gh']), 'flag'] = 2
df2.loc[(df2['flag']>-1) & (df2['flag']<1) & (df2['DIF']>df2['erl_dif']), 'flag'] = 2

** added 2022-11-08 **

# diafores vlakeis graphs akolouthoun
# plotting of GHI and SUMW values

# values in the same graph
# versus time (the index)
plt.scatter(df2.index[df2['Z']<80], df2['GH'][df2['Z']<80], s=0.005, c='k')
plt.scatter(df2.index[df2['Z']<80], df2['sumw'][df2['Z']<80], s=0.005, c='r')
plt.title('GHI & DNcosZ+DIF')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['GHI', 'DNcosZ+DIF'], markerscale=100, loc='upper right')
plt.show()
# versus the zenith angle
plt.scatter(df2['Z'][df2['Z']<80], df2['GH'][df2['Z']<80], s=0.002, c='k')
plt.scatter(df2['Z'][df2['Z']<80], df2['sumw'][df2['Z']<80], s=0.002, c='r')
plt.title('GHI & DNcosZ+DIF')
plt.ylabel('Irradiance [$W/m^2$]')
plt.xlabel('Zenith Angle [°]')
plt.legend(['GHI', 'DNcosZ+DIF'], markerscale=150, loc='upper right')
plt.show()

# as a difference (GHI-SUMW)
# versus time (the index)
plt.scatter(df2.index[df2['Z']<80],
            (df2['GH'][df2['Z']<80]-df2['sumw'][df2['Z']<80]), s=0.005, c='k')
plt.hlines(y=0, xmin=df2.index.min(),
           xmax=df2.index.max(), color='r')
plt.title('GHI-SUMW')
plt.ylabel('GHI-SUMW Difference [$W/m^2$]')
plt.show()
# versus the zenith angle
plt.scatter(df2['Z'][df2['Z']<80],
            (df2['GH'][df2['Z']<80]-df2['sumw'][df2['Z']<80]), s=0.005, c='k')
plt.hlines(y=0, xmin=10, xmax=85, color='r')
plt.title('GHI-SUMW')
plt.ylabel('GHI-SUMW Difference [$W/m^2$]')
plt.xlabel('Zenith Angle [°]')
plt.show()

plt.scatter(df2.index[df2['flag']==0], df2['GH'][df2['flag']==0], s=0.005, c='k')
plt.show()
'''
