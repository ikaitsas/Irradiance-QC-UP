# Irradiance-QC-UP
Quality control for the data of the radiometric station of the University of Patras.

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
This value is take from the 'Fysikh Atmosfairas II' course.
The formula for its variation through the year is:
Sa=Sav[1+0.033cos(360n/365)]   (from the 'Systhmata Hliakhs Energeias' course)

All negative GH(mV) or DIF(mV) values are set to 0. And rejected (flag==-1).
Datapoints that have Z>80 are also flagged with -1.
Datapoints with flag==-1 do not undergo any testing, as they are rejected.
Negative DN values are not flagged or corrected (this might need change).
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
ratio test, the diffuse ratio test, or both. Might change flag values to
differentiate betwwen the 2 tests.

The test results and limit values are all stored in dataframe "df2". It is then
saved in a .csv file, as it is. The file's/dataframe's column names mean:
UTC: datetime in UTC - corresponds to a datapoint
DIF: diffuse horizontal irradiance
GH: gobal horizontal irradiance
DN: direct normal irradiance
flag: the datapoint flag, the meaning of each value is mentioned above
Z: solar zenith angle
m0: cosine of Z (cosZ)
Day Count: Number of the day in the year - January 1st = 1
Sa: extraterrestrial irradiance (at the top of the atmosphere). Day Count is
    used for its calculation through the year
ppl_gh: physical possible limit for GH
ppl_dif: physical possible limit for DIF
erl_gh: extremely rare limit for GH
erl_dif: extremely rare limit for DIF
sumw: DN * cosZ + DIF, as calculated from the data
closr: closure ratio - the closure equation is GH = DNcosZ+DIF, with the ratio
       here being closr=GH/sumw, calculated from the data
dif_r: diffuse ratio - DIF/GH

In the .csv file, NaN values are stored as '', which are empty(?) spaces. That
might be changed to something like 'NaN' or 'null', I'll see. It needs to be
specified in what string/form Nan values are stored probably, maybe in the
_stats file?? I ideally want columns to not be objects...
