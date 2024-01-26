#Ashley Elliott
#DSEP isochrone data for HD 29391
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
#import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d

#%%
#Reading in each year
data20 = pd.read_csv("DSEP_20Myr.txt", sep = '\s+')
data25 = pd.read_csv("DSEP_25Myr.txt", sep = '\s+')
data30 = pd.read_csv("DSEP_30Myr.txt", sep = '\s+')
data315 = pd.read_csv("DSEP_31_5Myr.txt", sep = '\s+')
data32 = pd.read_csv("DSEP_32Myr.txt", sep = '\s+')
data325 = pd.read_csv("DSEP_32_5Myr.txt", sep = '\s+')
data33 = pd.read_csv("DSEP_33Myr.txt", sep = '\s+')
data335 = pd.read_csv("DSEP_33_5Myr.txt", sep = '\s+')
data34 = pd.read_csv("DSEP_34Myr.txt", sep = '\s+')
data345 = pd.read_csv("DSEP_34_5Myr.txt", sep = '\s+')
data35 = pd.read_csv("DSEP_35Myr.txt", sep = '\s+')

#%%
#reading in data we care about
log_L = []
log_T = []
L_20 = data20['LogL/Lo'].values.tolist()
log_L.append(L_20)
L_25 = data25['LogL/Lo'].values.tolist()
log_L.append(L_25)
L_30 = data30['LogL/Lo'].values.tolist()
log_L.append(L_30)
L_315 = data315['LogL/Lo'].values.tolist()
log_L.append(L_315)
L_32 = data32['LogL/Lo'].values.tolist()
log_L.append(L_32)
L_325 = data325['LogL/Lo'].values.tolist()
log_L.append(L_325)
L_33 = data33['LogL/Lo'].values.tolist()
log_L.append(L_33)
L_335 = data335['LogL/Lo'].values.tolist()
log_L.append(L_335)
L_34 = data34['LogL/Lo'].values.tolist()
log_L.append(L_34)
L_345 = data345['LogL/Lo'].values.tolist()
log_L.append(L_345)
L_35 = data35['LogL/Lo'].values.tolist()
log_L.append(L_35)

T_20 = data20['LogTeff'].values.tolist()
log_T.append(T_20)
T_25 = data25['LogTeff'].values.tolist()
log_T.append(T_25)
T_30 = data30['LogTeff'].values.tolist()
log_T.append(T_30)
T_315 = data315['LogTeff'].values.tolist()
log_T.append(T_315)
T_32 = data32['LogTeff'].values.tolist()
log_T.append(T_32)
T_325 = data325['LogTeff'].values.tolist()
log_T.append(T_325)
T_33 = data33['LogTeff'].values.tolist()
log_T.append(T_33)
T_335 = data335['LogTeff'].values.tolist()
log_T.append(T_335)
T_34 = data34['LogTeff'].values.tolist()
log_T.append(T_34)
T_345 = data345['LogTeff'].values.tolist()
log_T.append(T_345)
T_35 = data35['LogTeff'].values.tolist()
log_T.append(T_35)

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)
age = [20e7,25e7,30e7,31.5e7,32e7,32.5e7,33e7,33.5e7,34e7,34.5e7,35e7]

for x in range(0,11,1):
    #print('Luminosity for:',age[x])
    #print(log_L[x])
    plt.plot(log_T[x], log_L[x], '-.', label = f"{age[x]:.3e} years", markersize = 2, linewidth = 0.5)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    #plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()