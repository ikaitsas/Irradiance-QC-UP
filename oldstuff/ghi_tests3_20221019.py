# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

This code is an attempt to perform quality control in a solar irradiance
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


** added 2022/10/19 **

The Earth-Sun distance adjusted mean solar constant is Sav=1366 W/m^2.
H timh auth einai parmenh apo tis diafaneies tou 'Fysikh Atmosfairas II'.
The formula for its variation through the year is:
Sa=Sav[1+0.033cos(360n/365)]   (apo 'Systhmata Hliakhs Energeias')

Ta diagrammata tha ta peripoihthw kapoia stigmh sto mellon.

@author: yiann
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pytz


# latitude (lat) and longitude (lon)
lat = 38.29136389  # (>0 north & <0 south)
lon = 21.78858333  # (>0 east & <0 west)
# longitude of location's standard meridian (prime meridian at 0)
lon_m = 30.0  # Greece is +2 hours from 0 (15deg per standard meridian)


#%% Zenith Angle Calculation


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

# making the index a UTC+00 Datetime (Datetime is in local time - EET/UTC+02)  
df.index = df['Datetime']
# df = df.rename_axis(index = 'Datetime UTC')
# df.index.rename('Daytime UTC', inplace=True)
df.index.name = 'Datetime UTC'
grc_std = pytz.timezone('Etc/GMT-2')  # EET (UTC+02)
utc = pytz.timezone('UCT')  # UTC+00
df.index = df.index.tz_localize(grc_std).tz_convert(utc)



#%% GHI, DNI, DIF, Sa and Zenith Angle Values Input


# irradiance values input  
df2 = pd.read_csv(data_file, index_col=[0], usecols=[0, 4, 6, 8], sep=',', header=None,
                  parse_dates=True, na_values='"NAN"')
df2.columns = ['DIF', 'GH', 'DN']  # diffuse-global horizontal-direct normal
df2['DIF'] = df2['DIF'] * 1000 / 8.64  # Converts mV to W/m^2
df2['GH'] = df2['GH'] * 1000 / 8.63  # Converts mV to W/m^2
df2.index = df2.index.tz_localize('UTC')  # state that index is in UTC

df2 = df2.join(df['Z'])  # left join is default
df2 = df2.join(df['cosZ'])
df2['m0'] = df2['cosZ'].clip(lower=0)  # m0(=cosZ)>0 (or Z<90) for the tests
df2['BI'] = df2['DN'] * df2['m0']  # direct horizontal

# drop duplicates (meta ta join PANTA)
df2=df2[~df2.index.duplicated(keep='first')]  # or keep='last'??
'''
** added 2022/10/20 **

ta NaN sto telos mpainoun logw ths mikrhs anantistoixias??
sthn arxh, elegxontas to df.index.duplicated().sum(), exw 5 duplicates
ta duplicate indexes einai ola sto 2021-4-29 13:56-14:00 euros sto df
meta ta join ginontai 20 ta duplicated indexes (15+5) {indices to swsto..}
sto df2 exw sto 13:56-14:00 apo 3 duplicates kai sto 15:56-16:00 apo ena
kathe join prosthetei akoma mia timh?? alla meta giati ta akyra duplicates??
ti symvainei?? GIATI??
ta indexes akolouthoun thn seira tou arxeiou sto df, xronolikh seira sto df2
oi times ghi, dif, dni einai idies apo oso eida kai sta 5 duplicated indexes
se sxesh me ta arxika tous, opote den tithetai thema diotrhwshs autwn twn
timwn tou arxeiou, logika??

ta duplicated indexes exoun proelthei apo duplicated datetimes sto arxiko
arxeio (tha tsekarw an diorthenetai allh sthlh se auto kai se alla)
mporw me asfaleia(?) na na afairesw ta duplicate indexes ypo auto to context
me ena aplo df[~df.index.duplicated()] sta df, df2 eimai logika kalymmenos?
na valw keep='fisrt' omws h oxi?? tha rwthsw..

ta duplicated values genikotera vrhka pws dinontai apo:
count_series=df2.pivot_table(columns=['Column Name'], aggfunc='size')
count_series[count_series>1]
me index ths count_series thn sthlh me onoma 'Column Name'
tha yparxei kai pio apodotikos tropos, alla auton skefthka twra..

na frontisw na einai h metatroph anexarthth tou sygkekrimenou arxeiou kai na
pragmatopoieitai gia kathe eisagwmeno arxeio me ton idio tropo
an sta alla arxeia ta duplicate datetimes exoun shmasia prepei na to lavw
ypopsin mou kai auto logika
tha rwthsw ton k. argyriou, pros to paron to afhnw etsi, me ta duplicates..  
'''

# solar constant calculation (Sav = 1366 W/m2)
df2['Day Count'] = df2.index.dayofyear
'''
df2.index.size=525561 (gia 2021), enw exw 525600 lepta se mh disekto xrono
apo to df2['Day Count'].value_counts() vlepw apo pou leipoun metrhseis, afou
exw 24x60=1440 lepta mesa se mia hmera
p.x. leipoun 39 metrhseis sto 2021 arxeio
'''
df2['Sa'] = 1366 * ( 1 + 0.033 * np.cos( np.deg2rad((360 * df2['Day Count'])/365)) ) 
# alliws an thn thewrhsw statherh Sa=1366 W/m^2: df2['Sa'] = 1366 

'''
** added 2022-10-22 **

to: df2.drop(df2.loc[df2.index > '2021-12-31 21:59:00'].index, inplace=True)
diwxnei ta teleutaia datetimes pou einai NaN logw eisagwghs ths Z sto df2
to: df2.dropna() diwxnei ola ta pithana NaN stoixeia se kathe sthlh
sto prwto: size = 525441
sto deutero: size = 525423
akoma kai me to filtrarisma tou non-matching tail exw NaN values, sta cosZ, Z
df.isnull().values.any() na tsekarw an egine patata edw..
alla sto df den exw NaN values se auta.. ti symbainei??

vriskw ta datetime indices opou exw NaN times ws exhs:
index = df2['Z'].index[df2['Z'].apply(np.isnan)]
ti to prokalei twra auto??
to df.loc['2021-4-29 20:04:00'] den yparxei, alla to antistoixo df2.loc
yparxei, wtf??
mhpws kalytera na ftiakw ena date_range gia ton ypologismo twn Z, cosZ??

** added 2022-10-23 **

VRHKA TO PROBLHMAAAA, ORISTE H LYSHH????
to df2.loc['2021-4-29 22:04:00'] den yparxei, auto meiwnetai kata 2 wres sthn 
UTC metatroph tou df, opote kai xanetai h swsth antistoixia twn missing me
ta mh missing datetimes, mias kai to df.loc['2021-4-29 20:04:00'] fainetai
meta missing sta left joined df2, df['Z']..
apo edw paizei na proekypsan kai taparapanw duplicates??

*************************************************************************
* ARA NA YPOLOGISW TA cosZ, Z ME date_range() ENTOLH, GIA NA TO APOFYGW *
*************************************************************************
'''
df2 = df2.dropna()  # Deletes missing values (the tail moslty)



#%% Limits Tests


# Global Horizontal = GH
# DIffuse Horizontal = DIFF
# Direct Normal = DN
# Beam Irradiance = BI (BI=DNcosZ)

# PPL tests (QC1)
df2['ppl_min'] = -4  # minimun value for PPL test - common for all
df2['ppl_gh'] = df2['Sa'] * 1.5 * df2['m0']**1.2 + 100  
df2['ppl_dif'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 50  
df2['ppl_dn'] = df2['Sa']  
df2['ppl_bi'] = df2['Sa'] * df2['m0'] 

# ERL tests (QC2)
df2['erl_min'] = -2  # minimun value for ERL test - common for all
df2['erl_gh'] = df2['Sa'] * 1.2 * df2['m0']**1.2 + 50  
df2['erl_dif'] = df2['Sa'] * 0.75 * df2['m0']**1.2 + 30  
df2['erl_dn'] = df2['Sa'] * 0.95 * df2['m0']**0.2 + 10  
df2['erl_bi'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 10  # 

# plotting QC1 & QC2 (only values that correspond to Z<93 are shown)

# GH plots
plt.scatter(df2['Z'][df2['Z']<93], df2['ppl_gh'][df2['Z']<93], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<93], df2['erl_gh'][df2['Z']<93], s=0.001)
plt.scatter(df2['Z'][df2['Z']<93], df2['GH'][df2['Z']<93], s=0.001, c='k')
plt.show()

# DIF plots
plt.scatter(df2['Z'][df2['Z']<93], df2['ppl_dif'][df2['Z']<93], s=0.001, c='r') 
plt.scatter(df2['Z'][df2['Z']<93], df2['erl_dif'][df2['Z']<93], s=0.001)
plt.scatter(df2['Z'][df2['Z']<93], df2['DIF'][df2['Z']<93], s=0.001, c='k')
plt.show()

# DN plots
plt.scatter(df2['Z'][df2['Z']<93], df2['ppl_dn'][df2['Z']<93], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<93], df2['erl_dn'][df2['Z']<93], s=0.001)
plt.scatter(df2['Z'][df2['Z']<93], df2['DN'][df2['Z']<93], s=0.001, c='k')
plt.show()

# BI plots  # mporei na mhn xreiazontai??
plt.scatter(df2['Z'][df2['Z']<93], df2['ppl_bi'][df2['Z']<93], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<93], df2['erl_bi'][df2['Z']<93], s=0.001)
plt.scatter(df2['Z'][df2['Z']<93], df2['BI'][df2['Z']<93], s=0.001, c='k')
plt.show()



#%% Comparison Tests


# closure equation test
df2['sumw'] = df2['BI'] + df2['DIF']
df2['closr'] = df2['GH'] / df2['sumw']
# to ena paper grafei sumw/ghi, to allo (apo to bsrn) ghi/sumw??

# sto ena paper exei ghi>50 kai sto allo sumw>50, ti na valw? tha rwthsw..
plt.scatter( df2['Z'][ (df2['Z']<93) & (df2['sumw']>50) ],
            df2['closr'][ (df2['Z']<93) & (df2['sumw']>50) ], s=0.001, c='k')
'''
mallon na filtrarw apo prin me mask() h where() to closr ws pros to sumw,
gia ma einai pio euanagnwstos o kwdikas (for Z<93, sumw<50 test not possible)
apofeugw kai ta inf values me auto ton tropo
kai na kanw  kai to idio gia to diffuse ratio test
prepei na ginei giati thelw na kanw kai pass/fail pososta, opote kai na
exw kalytera mono ta katallhla datetime indexes gia auta
'''
# ratio limits visualization
plt.hlines(y=1.08, xmin=10, xmax=75, color='r')
plt.hlines(y=1.15, xmin=75, xmax=93, color='r')
plt.vlines(x=75, ymin=1.08, ymax=1.15, color='r')
plt.hlines(y=0.92, xmin=10, xmax=75, color='r')
plt.hlines(y=0.85, xmin=75, xmax=93, color='r')
plt.vlines(x=75, ymin=0.85, ymax=0.92, color='r')
plt.ylim([0,2])  # kapoia pane mexri to y=6..
plt.show()
'''
poly o,ti na nai vgainei to parapanw..
prepei na kanw kai ta pass/fail pososta na exw plhrh eikona..
mhpws epeidh einai pio 'mikra' ta DN/BI sta mikra Z apo to 'anamenomeno'..
san na fainetai na kanoun mia voutia sta mikra Z??
den paizei na peiraksa kati elpizw..
p.x. df2['DN'][df2['Z']<kapoia timh].median() goes up for Z going up
mallon vlakeies grafw..
'''
# diffuse ratio test
df2['dif_r'] = df2['DIF'] / df2['GH']

plt.scatter( df2['Z'][ (df2['Z']<93) & (df2['GH']>50) ],
            df2['dif_r'][ (df2['Z']<93) & (df2['GH']>50) ], s=0.001, c='k') 
# ratio limits visualization
plt.hlines(y=1.05, xmin=10, xmax=75, color='r')
plt.hlines(y=1.1, xmin=75, xmax=93, color='r')
plt.vlines(x=75, ymin=1.05, ymax=1.1, color='r')
plt.show()



# %% Test Plots

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
'''
# aplo plotarisma twn timwn aktinovolias 
plt.plot(df2['GH'])
plt.show()
plt.plot(df2['DN'])
plt.show()
plt.plot(df2['DIF'])
plt.show()
