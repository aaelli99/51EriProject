import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research")
import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d
import csv
import corner
import pandas as pd
# from isochrones.mist import MISTIsochroneGrid

#%%
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
iso = read_mist_models.ISO('MIST_iso_FEH013_0622_v2.iso')
#%%
#extracting the information from the MIST outputs
eepnum = []
age = []
log_Teff = []
log_L = []
log_r = []
star_mass = []
EEP_num = []

for i in range(iso.num_ages):
# for i in range(1):
    age.append(iso.ages[i])
    eepnum.append(iso.isos[iso.age_index(age[i])]['EEP'])
    
    for j in range(len(eepnum[i])):
        if eepnum[i][j] <= 454 and eepnum[i][j]>=202:
            EEP_num.append(eepnum[i][j])
            print("It worked.")
            print("EEP num:", eepnum[i][j])
    
    # print("EEP array:", EEP_num[i])

    lums = np.empty(len(EEP_num))
    rads = np.empty(len(EEP_num))
    temps = np.empty(len(EEP_num))
    mass = np.empty(len(EEP_num))
    for k in range(len(lums)):
        lums[k] = iso.isos[iso.age_index(age[i])]['log_L'][k]
        rads[k] = iso.isos[iso.age_index(age[i])]['log_R'][k]
        temps[k] = iso.isos[iso.age_index(age[i])]['log_Teff'][k]
        mass[k] = iso.isos[iso.age_index(age[i])]['star_mass'][k]
        
    log_L.append(lums)
    log_Teff.append(temps)
    log_r.append(rads)
    star_mass.append(mass)
    EEP_num = []
    
#%%
#Star information 
#new star info date: 06/27
L_star = np.log10(5.712)
Teff_star = np.log10(7424)
dL = 0.434*(0.096/5.712)
dT = 0.434*(45/7424)

#%%
for jj in range(0,80,1):
    plt.plot(log_Teff[jj], log_L[jj], '-s', label = f"{round(((age[jj]))/10**6,3)} Myr (MIST)", markersize = 3, linewidth = 0.5)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.5, 1])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 2, capsize = 3, color = 'k', label = 'HD29391')
# plt.legend(prop={'size': 15})
# for kk in range(24,34,1):
#     plt.plot(logTeff[kk], logL[kk], '-b', label = f"{round((10**(age_P[kk]))/10**6,3)} Myr (PARSEC)", markersize = 3, linewidth = 0.5)

plt.title('MIST Isochrones with [Fe/H] = 0.13')
# plt.legend(prop={'size': 15})
plt.show()

#%%
new_L = []
new_T = []
new_M = []
for ii in range(0,iso.num_ages,1):  

    f_age = interp1d(log_Teff[ii], log_L[ii])
    f_mass = interp1d(log_Teff[ii],star_mass[ii])
    new_temp = np.arange(min(log_Teff[ii]), max(log_Teff[ii]), 0.0001)
    new_lum = np.empty(len(new_temp))
    new_mass = np.empty(len(new_temp))
    new_lum = f_age(new_temp)
    new_mass = f_mass(new_temp)
    new_L.append(new_lum)
    new_T.append(new_temp)
    new_M.append(new_mass)
#%%
#Distance minimization code
def min_distance(x1,y1):

    min_dists = []
    min_idx = []
 
    for j in range(len(age)):   
        d = 10
        idx = 0
        dist = np.empty(len(new_T[j]))
        for i in range(len(new_T[j])):
            # d = []
            # idx = []
            x2 = new_T[j][i]
            y2 = new_L[j][i]
            #print("x2:", x2, "y2:", y2)
            dist[i] = np.sqrt((abs(x2-x1)**2)+(abs(y2-y1)**2))
            # print("Distance:", dist[i])
            #         #print("Old min dist:", d)
            if dist[i] < d:
                # print("Dist:",dist[i]," is less than ", d)
                d = dist[i]
                # print("Index position is:", i)
                idx = i
                # print("New minimum distance is:", d, "with index:", idx)
        min_dists.append(d)
        min_idx.append(idx)
            


    min_iso = np.min(min_dists)     #finds the shortest distance
    min_pos = np.argmin(min_dists)  #finds the position of the shortest distance
    min_age = age[min_pos]          #finds the age in the age array
    min_mass = new_M[min_pos][min_idx[min_pos]]       #finds the mass
        
    print("The closest isochrone is:", (min_age), "Myr") 
    # print("The minimum distance is:", min_iso) 
    # print("The index is:", min_pos)
    print("The mass is:", min_mass, "solar masses")
    
    return((min_age,min_mass))

#%%
#Monte Carlo simulation
num_iter = 1000
age_MC = []
mass_MC = []
newvals = []

for b in range(num_iter):
    
    L_err = np.random.normal(-dL,dL)
    T_err = np.random.normal(-dT,dT)
    new_Teff_star = Teff_star+T_err
    new_L_star = L_star+L_err
    #print("L:\n", L_array)
    #print("T:\n", T_array)
    #print("New L: ", new_L_star, "New T:", new_Teff_star)
    dist_calc = min_distance(new_Teff_star, new_L_star)
    
    newvals.append((new_Teff_star,new_L_star))    
    
    age_MC.append(dist_calc[0])
    mass_MC.append(dist_calc[1])
    
    print("Iteration ", b+1, "out of ", num_iter)
    
#%%
#converts data out of log space
new_temps = []
age_yr = []
lum_nl = []
for i in range(num_iter):
    t = 10**(newvals[i][0])
    a = ((age_MC[i]))/(1000000)
    l = 10**(newvals[i][1])
    new_temps.append(t)
    age_yr.append(a)
    lum_nl.append(l)
#%%
#plots the corner plots
samples = np.vstack([age_yr, mass_MC, new_temps, lum_nl])
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),range = [(15,60),(1.545,1.595), (7250,7500),(5.3,5.9)],labels = [r"$Age [Myr]$", r"$Mass [M_{\odot}]$", r"$T_{eff} [K]$", r"$L_{\star} [L_{\odot}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((4,4))

#%%
from numpy import savetxt
MC_L = np.array(lum_nl)
MC_T = np.array(new_temps)
MC_age = np.array(age_yr)
MC_mass = np.array(mass_MC)

dat_MC = np.array([MC_L, MC_T, MC_age, MC_mass])
dat_MC = dat_MC.T
savetxt('MIST_FEH013_0622_v2.csv', dat_MC, delimiter = ',')