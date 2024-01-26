#Ashley Elliott
#isochrone data for HD 29391
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d
import csv
import corner
#%%
# os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data\\MIST_iso_6331c71466188.iso")
#reading in the file
# iso = read_mist_models.ISO('MIST_iso_6331c71466188.iso')
iso = read_mist_models.ISO('MIST_iso_64920b005c8e5.iso')
#iso = read_mist_models.ISO('tmp1670874098.iso')
#%%
#checking header info
print('version: ', iso.version)
print('abundances: ', iso.abun)
print('rotation: ', iso.rot)
print('ages: ', [round(x,2) for x in iso.ages])
print('number of ages: ', iso.num_ages)
print('available columns: ', iso.hdr_list)

#%%
#extracting the information from the MIST outputs
eepnum = []
age = []
log_Teff = []
log_L = []
log_r = []
star_mass = []
EEP_num = []

for i in range(41):
    age.append(iso.ages[i])
    eepnum.append(iso.isos[iso.age_index(age[i])]['EEP'])
    
    for j in range(len(eepnum[i])):
        if eepnum[i][j] <= 300:
            EEP_num.append(eepnum[i][j])

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
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)

# L_star = 5.712
# Teff_star = 7414
# dL = 0.096
# dT = 32
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


for x in range(45,55,1):
    plt.plot(log_Teff[x], log_L[x], '-s', label = f"{10**(age[x]):.2e} years", markersize = 5, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%

x1 = L_star
y1 = Teff_star
dist = []
min_dist = []
for m in range(0,41,1):
    chrone_L = log_L[m]
    chrone_T = log_Teff[m]
    d = np.empty(len(chrone_L))
    for k in range(len(chrone_L)):
        x2 = chrone_L[k]
        y2 = chrone_T[k]
        d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
    
    dist.append(d)
    
for n in range(len(dist)):
    min_dist.append(np.min(dist[n]))

min_iso = np.min(min_dist)
pos_iso = np.argmin(min_dist)

print("The minimum distance between points is", min_iso)
print("Age of isochrone:",f"{10**(age[pos_iso]):.2e} years")

    

#%%
plt.plot(log_Teff[pos_iso], log_L[pos_iso], '-s', label = f"{10**(age[pos_iso]):.2e} years", markersize = 1.0, linewidth = 0.2)
plt.xlabel('log(Teff)')
plt.ylabel('log(L)')
plt.axis([4,3.7,0,1.5])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%
Teff = []
Lum = []
for v in range(0,107,1):
    tem = 10**log_Teff[v]
    bright = 10**log_L[v]
    
    Lum.append(bright)
    Teff.append(tem)

        

#%%
#Interpolation using temperature as the independent variable
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


new_L = []
new_T = []
new_M = []
for ii in range(0,41,1):  

    f_age = interp1d(log_Teff[ii], log_L[ii])
    f_mass = interp1d(log_Teff[ii],star_mass[ii])
    new_temp = np.arange(min(log_Teff[ii]), max(log_Teff[ii]), 0.0001)
    # f_age = interp1d(Teff[ii], Lum[ii])
    # f_mass = interp1d(Teff[ii],star_mass[ii])
    # new_temp = np.arange(min(Teff[ii]), max(Teff[ii]), 0.1)
    new_lum = np.empty(len(new_temp))
    new_mass = np.empty(len(new_temp))
    new_lum = f_age(new_temp)
    new_mass = f_mass(new_temp)
    new_L.append(new_lum)
    new_T.append(new_temp)
    new_M.append(new_mass)
    
#%%
#Interpolation using mass as the independent variable instead of temperature
#April 18,2023
new_L = []
new_T = []
new_M = []
for ii in range(0,107,1):  

    f_massL = interp1d(star_mass[ii], log_L[ii])
    f_massT = interp1d(star_mass[ii], log_Teff[ii])
    
    new_mass = np.arange(min(star_mass[ii]), max(star_mass[ii]),0.001)
    new_lum = np.empty(len(new_mass))
    new_temp = np.empty(len(new_mass))
    new_lum = f_massL(new_mass)
    new_temp = f_massT(new_mass)
    new_L.append(new_lum)
    new_T.append(new_temp)
    new_M.append(new_mass)



#%%
for jj in range(45,55,1):
    plt.plot(new_T[jj], new_L[jj], '-s', label = f"{round((10**(age[jj]))/10**6,1)} Myr", markersize = 3, linewidth = 0.5)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 2, capsize = 3, color = 'k', label = 'HD29391')
plt.legend(prop={'size': 15})
plt.title('MIST Isochrones')
plt.show()

#%%
y1 = L_star
x1 = Teff_star
dist = []
min_dist = []
for m in range(0,107,1):
    chrone_L = new_L[m]
    chrone_T = new_T[m]
    d = np.empty(len(chrone_L))
    for k in range(len(chrone_L)):
        y2 = chrone_L[k]
        x2 = chrone_T[k]
        d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
    
    dist.append(d)
    
for n in range(len(dist)):
    min_dist.append(np.min(dist[n]))

min_iso = np.min(min_dist)
pos_iso = np.argmin(min_dist)

L = new_L[pos_iso]
T = new_T[pos_iso]
M = new_M[pos_iso]

l = np.asarray(L)
t = np.asarray(T)
idx_l = (np.abs(l-L_star)).argmin()
idx_t = (np.abs(t-Teff_star)).argmin()

print('Closest luminosity to ', L_star, ':', l[idx_l], 'position:', idx_l)
print('Closest temperature to ', Teff_star, ':', t[idx_t], 'position:', idx_t)


print("The minimum distance between points is", min_iso)
print("Age of isochrone:",f"{10**(age[pos_iso]):.2e} years")
print("Mass of the star (closest luminosity):", M[idx_l], 'Solar masses')
print("Mass of the star (closest temperature):", M[idx_t], 'Solar masses')

#%%
#minimum distance function

def min_distance(x1,y1,T_array,L_array):
    dist = []
    min_dist = []
    # chrone_L = L_array
    # chrone_T = T_array
    # d = np.empty(len(chrone_L))
    # for k in range(len(chrone_L)):
    #     x2 = chrone_T[k]
    #     y2 = chrone_L[k]
    #     d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
    # dist.append(d)
                     
    for m in range(0,41,1):
        chrone_L = L_array[m]
        chrone_T = T_array[m]
        d = np.empty(len(chrone_L))
        for k in range(len(chrone_L)):
            y2 = chrone_L[k]
            x2 = chrone_T[k]
            d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
    
        dist.append(d)
    
    for n in range(len(dist)):
        min_dist.append(np.min(dist[n]))

    min_iso = np.min(min_dist)
    pos_iso = np.argmin(min_dist)

    L = new_L[pos_iso]
    T = new_T[pos_iso]
    M = new_M[pos_iso]
    
    l = np.asarray(L)
    t = np.asarray(T)
    idx_l = (np.abs(l-L_star)).argmin()
    idx_t = (np.abs(t-Teff_star)).argmin()
    
    age_iso = age[pos_iso]
    mass_iso_L = M[idx_l]
    mass_iso_T = M[idx_t]
    
    return((age_iso, mass_iso_L, mass_iso_T, L, T, M))
#%%
#Monte Carlo simulations
num_iter = 500
best_L = []
best_T = []
best_M = []
age_MC = []
massL_MC = []
massT_MC = []
newvals = []

for i in range(num_iter):
    
    L_err = np.random.normal(-dL,dL)
    T_err = np.random.normal(-dT,dT)
    new_Teff_star = Teff_star+T_err
    new_L_star = L_star+L_err
    
    dist_calc = min_distance(new_Teff_star, new_L_star, new_T, new_L)
    
    newvals.append((new_Teff_star,new_L_star))    
    
    best_L.append(dist_calc[3])
    best_T.append(dist_calc[4])
    best_M.append(dist_calc[5])
    age_MC.append(dist_calc[0])
    massL_MC.append(dist_calc[1])
    massT_MC.append(dist_calc[2])
    
    print("Iteration ", i+1, "out of ", num_iter)

#%%
for i in range(num_iter):
    plt.plot(newvals[i][0],newvals[i][1],'.')
#%%
list_columns = ["MC Age", "MC Mass_L", "MC Mass_T", "T_star val", "L_star val"]

file = open("MIST_MC_v12.csv", "w", newline ='')
writer = csv.writer(file)
writer.writerow(list_columns)
for w in range(len(age_MC)):
    writer.writerow([age_MC[w], massL_MC[w], massT_MC[w], newvals[w][0], newvals[w][1]])
    
file.close()


#%%
new_temps = []
age_yr = []
lum_nl = []
for i in range(num_iter):
    t = 10**(newvals[i][0])
    a = (10**(age_MC[i]))/(1e6)
    l = 10**(newvals[i][1])
    new_temps.append(t)
    age_yr.append(a)
    lum_nl.append(l)
#%%
samples = np.vstack([age_yr, massT_MC, new_temps, lum_nl])

mean_age = np.mean(age_MC)
mean_mass = np.mean(massT_MC)
std_age = np.std(age_MC)
std_mass = np.std(massT_MC)

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


figure = corner.corner(samples.T,quantiles=(0.15, 0.5,0.85),labels = [r"$Age [Myr]$", r"$Mass [M_{\odot}]$", r"$T_{eff} [K]$", r"$L_{\star} [L_{\odot}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((4,4))

value1 = np.mean(samples.T,axis=0)

# for i in range(4):
#     ax = axes[i,i]
#     ax.axvline(value1[i], color = "r")
    

# for yi in range(4):
#     for xi in range(yi):
#         ax1 = axes[yi, xi]
#         ax1.axvline(value1[xi], color="r")
#         ax1.axhline(value1[yi], color="r")
#         ax1.plot(value1[xi], value1[yi], "sr")
 
# ax2 = axes[1,1]
# ax2.set_xlim()