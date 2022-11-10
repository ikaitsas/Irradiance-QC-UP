# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:28:42 2022

This code is an attempt to perform quality control οn a solar irradiance
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

# import irradiance data and create datetimes for zenith angle calculation
data_file = 'Solar_1min_2021.txt'
dtr=pd.date_range(start='2021-1-1 00:00:00', end='2022-1-1 1:59:00',
                 freq='min')



#%% Zenith Angle Calculation

# create datetime dataframe for specified year
df=pd.DataFrame(dtr, columns=['Datetime'])

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

# Global Horizontal = GH
# DIffuse Horizontal = DIF
# Direct Normal = DN
# Beam Irradiance = BI (BI=DNcosZ)

# irradiance values input
df2 = pd.read_csv(data_file, index_col=[0], usecols=[0, 4, 6, 8], sep=',',
                  header=None, parse_dates=True, na_values='"NAN"')
df2.columns = ['DIF', 'GH', 'DN']
df2.index = df2.index.tz_localize('UTC')  # state that index is in UTC
df2['GH'] = df2['GH'].clip(lower=0)  # remove negative voltage values
df2['DIF'] = df2['DIF'].clip(lower=0)
# mporei na xreiazontai clip kai oi arnhtikes DN
df2['DIF'] = df2['DIF'] * 1000 / 8.64  # Converts mV to W/m^2
df2['GH'] = df2['GH'] * 1000 / 8.63

df2 = df2.join(df['Z'])  # left join is default
#df2 = df2[df2['Z']<80]  # keep Z<80 values
df2['m0'] = np.cos(np.deg2rad(df2['Z'])).clip(lower=0)  # m0=cosZ>0

# drop duplicate rows (meta ta join kalytera)
df2=df2[~df2.index.duplicated(keep='first')]

# solar constant calculation (Sav = 1366 W/m2)
df2['Day Count'] = df2.index.dayofyear
df2['Sa'] = 1366 * ( 1 + 0.033 * np.cos( np.deg2rad((360 * df2['Day Count'])/365)) )
# alliws an thn thewrhsw statherh (Sa=1366 W/m^2): df2['Sa'] = 1366

#df2 = df2.dropna()  # Deletes NaN/missing values



#%% Limit Tests

# PPL tests (QC1)
df2['ppl_gh'] = df2['Sa'] * 1.5 * df2['m0']**1.2 + 100
df2['ppl_dif'] = df2['Sa'] * 0.95 * df2['m0']**1.2 + 50

# ERL tests (QC2)
df2['erl_gh'] = df2['Sa'] * 1.2 * df2['m0']**1.2 + 50
df2['erl_dif'] = df2['Sa'] * 0.75 * df2['m0']**1.2 + 30

# plotting QC1 & QC2 (only values that correspond to Z<80 are shown)

# GH plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_gh'][df2['Z']<80], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_gh'][df2['Z']<80], s=0.001)
plt.scatter(df2['Z'][df2['Z']<80], df2['GH'][df2['Z']<80], s=0.001, c='k')
plt.title('Global Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Dataset'], markerscale=150)
plt.show()

# DIF plots
plt.scatter(df2['Z'][df2['Z']<80], df2['ppl_dif'][df2['Z']<80], s=0.001, c='r')
plt.scatter(df2['Z'][df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.001)
plt.scatter(df2['Z'][df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.001, c='k')
plt.title('Diffuse Horizontal Irradiance')
plt.xlabel('Zenith Angle [°]')
plt.ylabel('Irradiance [$W/m^2$]')
plt.legend(['PPL', 'ERL', 'Dataset'], markerscale=150)
plt.show()



#%% Comparison Tests

# closure equation test
df2['sumw'] = df2['DN'] * df2['m0'] + df2['DIF']
df2['closr'] = df2['GH'] / df2['sumw']

# sto ena paper exei ghi>50 kai sto allo sumw>50, ti na valw? tha rwthsw..
plt.scatter(df2['Z'][(df2['Z']<80) & (df2['sumw']>50)],
            df2['closr'][(df2['Z']<80) & (df2['sumw']>50)], s=0.001, c='k')
'''
mallon na filtrarw apo prin me mask() h where() to closr ws pros to sumw,
gia ma einai pio euanagnwstos o kwdikas (for Z<93, sumw<50 test not possible)
apofeugw kai ta inf values (pou ektos twn eligible values exarxhs vevaia)
kai na kanw  kai to idio gia to diffuse ratio test
prepei na ginei giati thelw na kanw kai pass/fail pososta, opote kai na
exw kalytera mono ta katallhla datetime indexes gia auta??
'''
# ratio limits visualization
plt.hlines(y=1.08, xmin=10, xmax=75, color='r')
plt.hlines(y=1.15, xmin=75, xmax=85, color='r')
plt.vlines(x=75, ymin=1.08, ymax=1.15, color='r')
plt.hlines(y=0.92, xmin=10, xmax=75, color='r')
plt.hlines(y=0.85, xmin=75, xmax=85, color='r')
plt.vlines(x=75, ymin=0.85, ymax=0.92, color='r')
plt.ylim([0,2])  # kapoia pane mexri to y=6..
plt.title('Closure Ratio (GH/SUMW)')
plt.xlabel('Zenith Angle [°]')
plt.show()
'''
poly o,ti na nai vgainei to parapanw..
prepei na kanw kai ta pass/fail pososta na exw plhrh eikona..
mhpws epeidh einai pio 'mikra' ta DN/BI sta mikra Z apo to 'anamenomeno'..
san na fainetai na kanoun mia voutia sta mikra Z??
den paizei na peiraksa kati elpizw..
p.x. df2['DN'][df2['Z']<kapoia timh].median() goes up for Z going up
mallon vlakeies grafw..
to oti einai estimate paizei kapoio rolo??
'''

# diffuse ratio test
df2['dif_r'] = df2['DIF'] / df2['GH']

plt.scatter(df2['Z'][(df2['Z']<80) & (df2['GH']>50)],
            df2['dif_r'][(df2['Z']<80) & (df2['GH']>50)], s=0.001, c='k')
# ratio limits visualization
plt.hlines(y=1.05, xmin=10, xmax=75, color='r')
plt.hlines(y=1.1, xmin=75, xmax=85, color='r')
plt.vlines(x=75, ymin=1.05, ymax=1.1, color='r')
plt.title('Diffuse Ratio (DIF/GH)')
plt.xlabel('Zenith Angle [°]')
plt.show()



#%% ti menei - sxolia

'''
na apothikeusw ta irradiance values (kai ta Z??), mazi me ta flags pou
prokyptoyn apo tous parapanw elegxous (p.x. pass, fail, suspect etc..)

na kanw sort to datetime index to df2 prwta kalou kakou prin ta join

na eksagw statistika gia ta pososta epityxias/apotyxias
na to kanw ws pososto olwn twn timwn h twn eligible gia kathe test mono??
tha rwthsw..

** added 2022/10/20 **

ta NaN sto telos mpainoun logw ths mikrhs anantistoixias??
an den kanw drop ta duplicate indices, an yparxoun:
elegxontas to df2.index.duplicated().sum(), exw 5 duplicates,
sto 15:56-16:00, gia to 2021-4-29
fainetai na eginan duplicate ta Z, cosZ (??)
poly akyra mphkan auta.. apo pou?? ti symvainei?? giati??

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
'''

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

plt.scatter(df2.index[df2['Z']<80], df2['ppl_gh'][df2['Z']<80],
            s=0.0002, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_gh'][df2['Z']<80], s=0.0002)
plt.scatter(df2.index[df2['Z']<80], df2['GH'][df2['Z']<80], s=0.005, c='k')
plt.title('GHI')
plt.ylabel('Irradiance [$W/m^2$]')
plt.show()

plt.scatter(df2.index[df2['Z']<80], df2['DN'][df2['Z']<80], s=0.005, c='k')
plt.title('DNI')
plt.ylabel('Irradiance [$W/m^2$]')
plt.show()

plt.scatter(df2.index[df2['Z']<80], df2['ppl_dif'][df2['Z']<80],
            s=0.0002, c='r')
plt.scatter(df2.index[df2['Z']<80], df2['erl_dif'][df2['Z']<80], s=0.0002)
plt.scatter(df2.index[df2['Z']<80], df2['DIF'][df2['Z']<80], s=0.005, c='k')
plt.title('DIF')
plt.ylabel('Irradiance [$W/m^2$]')
plt.show()

'''
idea: na kanw 2d diagramma timwn, p.x. zenitheia-hmeromhnia me xrwma gia
      thn optikopoihsh ths entashs ths aktinobolias
sta parapanw den einai eudiakrito to an ta data paraviazoun ta oria, logw tou
megalou plithous timwn, ena 2d diagramma isws tha voithouse se mia me to
mati epithewrhsh twn timwn pisteuw..
gia plaka kai eksoikeiwsh, oxi sta plaisia ths ergasias..
'''
