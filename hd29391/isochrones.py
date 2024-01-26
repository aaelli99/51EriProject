#Ashley Elliott
#isochrone data for HD 29391


import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d

#%%
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data\\MIST_iso_6331c71466188.iso")
#reading in the file
iso = read_mist_models.ISO('MIST_iso_6331c71466188.iso')

#%%
#checking header info
print('version: ', iso.version)
print('abundances: ', iso.abun)
print('rotation: ', iso.rot)
print('ages: ', [round(x,2) for x in iso.ages])
print('number of ages: ', iso.num_ages)
print('available columns: ', iso.hdr_list)


#%%
age = []
log_Teff = []
log_L = []
EEP = []
log_r = []
star_mass = []
EEP_Val = []

good = []
bad = []
gcounter = 0
bcounter = 0
for i in range(iso.num_ages):
    age.append(iso.ages[i])
    EEP.append(iso.isos[iso.age_index(age[i])]['EEP'])
    #removing the pre-main sequence star data
    for j in range(len(EEP[i])):
        if EEP[i][j] > 202:
            print("It worked")
            log_Teff.append(iso.isos[iso.age_index(age[i])]['log_Teff'][j:])
            log_L.append(iso.isos[iso.age_index(age[i])]['log_L'][j:])
            log_r.append(iso.isos[iso.age_index(age[i])]['log_R'][j:])
            star_mass.append(iso.isos[iso.age_index(age[i])]['star_mass'][j:])
            EEP_Val.append(EEP[i][j])
            gcounter+=1
            #break 
        else:
            bcounter+=1
            print("Too big")
#%%
age_list = []
feh_list = []
for x in range(len(star_mass)):
    age_list.append(np.ones(len(star_mass[x]))*age[x])
    feh_list.append(np.ones(len(star_mass[x]))*0.1)




# age_listf = age_list.flatten()
# feh_listf = feh_list.flatten()
# star_massf = star_mass.flatten()
# log_Tefff = log_Teff.flatten()
# log_Lf = log_L.flatten()
# log_rf = log_r.flatten()

age_listf = np.concatenate(age_list).ravel()
feh_listf = np.concatenate(feh_list).ravel()
star_massf = np.concatenate(star_mass).ravel()
log_Tefff = np.concatenate(log_Teff).ravel()
log_Lf = np.concatenate(log_L).ravel()
log_rf = np.concatenate(log_r).ravel()
    
inputs = list(zip(age_listf, star_massf, feh_listf))
outputs= list(zip(log_Tefff, log_Lf, log_rf))

f = LinearNDInterpolator(inputs,outputs)

#%%
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)

for x in range(0,25,5):
    plt.plot(log_Teff[x], log_L[x], '-s', label = f"{10**(age[x]):.2e} years", markersize = 1.0, linewidth = 0.2)
    plt.xlabel('log(Teff)')
    plt.ylabel('log(L)')
    #plt.axis([4.5, 3.5, 0.15, 2.25])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%
plt.plot(log_Teff[45], log_L[45], label = f"{10**(age[x]):.2e}")
plt.plot(Teff_star, L_star, "*")
plt.xlabel('log(Teff)')
plt.ylabel('log(L)')
plt.axis([6, 3.3, -5, 7])
plt.legend()
plt.show()

#%%
