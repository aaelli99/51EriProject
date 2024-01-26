#Planet mass evolutions

#%%
#Importing all necessary packages
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Planet")
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy import interpolate
import pandas as pd
from scipy.interpolate import interp1d
import random
import corner
import re
from scipy.interpolate import Rbf
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp2d
from scipy.interpolate import NearestNDInterpolator
#%%
#Reading in the data file and rewriting it to remove anything unneccessary
#only to be ran once on any new evolution table
f = open('mass_evolution.txt','r')
f_new = open('mass_evo_edit.txt','a')
# f_new.write("M/Msun  age(Gyr) log L/Lsun  Teff(K)  log g  R/Rsun   log I\n")
no = "#"
for line in f:
    if not no in line:
        f_new.write(line)
        f_new.write("\n")

f_new.close()      
f.close()

#%%
#Reading in the data into a pandas, sorting the data according to the masses to create individual curves
data = pd.read_csv("mass_evo_edit.txt", sep='\s+')

Age = (data['age(Gyr)']*(1e9))/(1e6)
mass = []
age = []
logL = []
T = []

m = [0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.0115,0.012,0.0125,0.013,0.0135,0.014,0.0145,0.015,0.0155,0.016,0.0165,0.017,0.0175,0.018,0.0185,0.019,0.0195,0.02,0.022,0.024,0.026,0.028,0.03,0.033,0.035,0.038,0.04,0.043,0.045,0.048,0.05,0.053,0.055,0.058,0.06,0.063,0.064,0.065,0.066,0.067,0.068,0.069,0.07,0.071,0.072,0.073,0.074,0.075,0.077,0.078,0.08]
lens = np.empty(len(m))

for ii in range(len(m)):
    idx = 0
    for i in range(len(data['M/Msun'])):
        if data["M/Msun"][i] == m[ii]:
            # print("Data is: ", data["logAge"][i], "with index ", i)
            # print("Age is: ", age[ii])
            idx = idx+1
            #print(idx)
        lens[ii] = idx
    
    beg_pos = np.empty(len(m))
    end_pos = np.empty(len(m))
    # print(lens[ii])
    tot = lens[0]-1
    beg_pos[0] = 0
    end_pos[0] = int(tot)
    for n in range(len(lens)-1):
        # print(x)
        beg_pos[n+1] = int(tot+1)
        tot = tot + lens[n+1]
        end_pos[n+1] = int(tot)
        
for nn in range(len(m)):
    #     print("Index:", y)
    #     print("Begin: ",beg_pos[y])
    #     print("End: ", end_pos[y])

    lums = data["logL/Lsun"][int(beg_pos[nn]):int(end_pos[nn])+1].values.tolist()
    temps = data["Teff(K)"][int(beg_pos[nn]):int(end_pos[nn])+1].values.tolist()
    a = Age[int(beg_pos[nn]):int(end_pos[nn])+1].values.tolist()
    pmass = data["M/Msun"][int(beg_pos[nn]):int(end_pos[nn])+1].values.tolist()  
    logL.append(lums)
    T.append(temps)
    age.append(a)
    mass.append(pmass)
    
#%%
#Constants needed for the plots and interpolation, includes age of planet with uncertainties, temperature of planet with uncertainties
#and turning the masses into jupiter masses instead of solar masses
age_planet = 23
da_up = 2.5
da_down = 1.5
age_error = [da_down, da_up]
# teff_planet = np.log10(760)
# dT = 0.434*(20/760)
teff_planet = 760
dT = 20
mass_J = []
for i in range(len(m)):
    M = (m[i]*(1.989e30))/(1.899e27)
    mass_J.append(M)
#%%
#Plotting the base mass evolution tracks
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

        
for r in range(0,10,1):
    plt.plot(age[r], T[r], '-s', label = f"{mass_J[r]:.2}" r'$M_{Jup}$', markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$Age [Myr]$', fontsize = 30)
    plt.ylabel(r'$Teff [K]$', fontsize = 30)
    plt.axis([0, 50, 0, 1000])
#plt.axis([3.9, 3.8,0.5,1.0])
plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.5], [2.5]]), yerr = dT, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = '51 Eri b')
plt.legend()
plt.show()

#%%
#Interpolates the mass tracks to get a finer grid spacing for the lines
new_a = []
new_T = []

for k in range(0,len(m),1):
    mass_func = interp1d(age[k],T[k], kind = 'cubic')
    new_age = np.linspace(min(age[k]), max(age[k]), 20000)
    new_temps = mass_func(new_age)
    new_T.append(new_temps)
    new_a.append(new_age)



#%%
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for it in range(0, 15,1):
    lab = (m[it]*(1.989e30))/(1.899e27)
    plt.plot(new_a[it], new_T[it], label = f"{lab:.3}" r'$M_{Jup}$', markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$Age [Myr]$', fontsize = 30)
    plt.ylabel(r'$T_{eff} [K]$', fontsize = 30)
    plt.axis([0,50,0,1000])
plt.errorbar(23,760, yerr = 20, xerr = np.array([[1.5],[2.5]]), linestyle = None, linewidth = 0.5, capsize = 3, color = 'black')
plt.legend()
plt.show()


#%%
#plotting the interpolated mass tracks
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


for q in range(0,15,1):
    plt.plot(new_a[q], new_T[q], '-s', label = f"{mass_J[q]:.3}" r'$M_{Jup}$', markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$Age [Myr]$', fontsize = 30)
    plt.ylabel(r'$Teff [K]$', fontsize = 30)
    plt.axis([0, 50, 0, 1000])

plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.5], [2.5]]), yerr = dT, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = '51 Eri b')
# plt.grid(which = 'major')
plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
plt.grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.legend()
plt.show()

#%%
#using an interpolator to interpolate between the given curves
# data = pd.read_csv("mass_evo_edit.txt", sep='\s+')
temp_p = data['Teff(K)']
age_p = data['age(Gyr)']
mass_p = data['M/Msun']

x = np.array(age_p)                             #age is in Gyr here
y = np.array(temp_p)                            #temperature is in Kelvin
z = np.array(mass_p)                            #masses are in solar masses here
age_planet_MC = 0.023                           #age is in Gyr
da_up_MC = 0.0025
da_down_MC = 0.0015
f = LinearNDInterpolator(list(zip(x,y)), z, rescale = True)


#%%
temp_p = data['Teff(K)']
age_p = data['age(Gyr)']
mass_p = data['M/Msun']

x = np.array(age_p)                             #age is in Gyr here
y = np.array(temp_p)                            #temperature is in Kelvin
z = np.array(mass_p)                            #masses are in solar masses here
# x = new_a.flatten()
# y = new_T.flatten()

new_f = interp2d(x,y, z, 'cubic')
age_planet_MC = 0.023                           #age is in Gyr
da_up_MC = 0.0025
da_down_MC = 0.0015

#%%
temp_p = data['Teff(K)']
age_p = data['age(Gyr)']
mass_p = data['M/Msun']

x = np.array(age_p)                             #age is in Gyr here
y = np.array(temp_p)                            #temperature is in Kelvin
z = np.array(mass_p)                            #masses are in solar masses here

interp = Rbf(x,y,z)

#%%
new_f(0.023,760)
#%%
mass_check = []

for g in range(len(x)):
    m_check = f(x[g],y[g])
    mc = np.around(m_check, 5)
    if mc== z[g]:
        print("It worked")
    else:
        print("Sad times.")
    mass_check.append(mc)
#%%

def assym_dist(neg,pos, val):
    r = np.random.normal(0,1)
    # print(r)
    if r<0:
        err = r*neg
        new_val = val+err
        # print("Value is negative.")
        # print("New value:", val, "+", err, "=", new_val)
    else:
        err = r*pos
        new_val = val+err
        # print("Value is positive.")
        # print("New value:", val, "+", err, "=", new_val)
    return(new_val)

#%%
#Monte Carlo simulation randomizing the uncertainties and then plugging into the interpolated function to output a 
#new mass

num_iter = 5000
mass_MC = []
newa = []
newt = []

for b in range(num_iter):
    A_err = assym_dist(da_down_MC, da_up_MC, age_planet_MC)
    T_err = np.random.normal(teff_planet,dT)
    # print("T_err:", T_err)
    # new_T_p = teff_planet+T_err
    # new_a_p = age_planet_MC+A_err

    new_p_mass = f(A_err, T_err)
    # newpmass = new_p_mass.tolist()
    print("Iteration ", b+1, "out of ", num_iter)
    
    newa.append((A_err*(1e9))/(1e6))
    newt.append(T_err)  
    mass_MC.append((new_p_mass*(1.989e30))/(1.899e27))
 

#%%
#plotting the corner plot for the Monte Carlo simulation

samples = np.vstack([newa, newt, mass_MC])
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 22})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),labels = [r"$Age [Myr]$", r"$T_{eff} [K]$",r"$Mass [M_{J}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((3,3))

#%%
# #%%
# #Distance minimization code
# def min_distance(x1,y1):

#     min_dists = []
#     min_idx = []
 
#     for j in range(len(age)):   
#         d = 10
#         idx = 0
#         dist = np.empty(len(new_T[j]))
#         for i in range(len(new_T[j])):
#             # d = []
#             # idx = []
#             x2 = new_a[j][i]
#             y2 = new_T[j][i]
#             #print("x2:", x2, "y2:", y2)
#             dist[i] = np.sqrt((abs(x2-x1)**2)+(abs(y2-y1)**2))
#             # print("Distance:", dist[i])
#             #         #print("Old min dist:", d)
#             if dist[i] < d:
#                 # print("Dist:",dist[i]," is less than ", d)
#                 d = dist[i]
#                 # print("Index position is:", i)
#                 idx = i
#                 # print("New minimum distance is:", d, "with index:", idx)
#         min_dists.append(d)
#         min_idx.append(idx)
            


#     min_iso = np.min(min_dists)     #finds the shortest distance
#     min_pos = np.argmin(min_dists)  #finds the position of the shortest distance
#     min_mass = mass[min_pos]          #finds the age in the age array
#     # min_mass = new_M[min_pos][min_idx[min_pos]]       #finds the mass
        
#     print("The closest mass is:",min_mass, "Jupiter masses") 
#     # print("The minimum distance is:", min_iso) 
#     # print("The index is:", min_pos)
#     # print("The mass is:", min_mass, "solar masses")
    
#     return(min_mass)

#%%
# #Monte Carlo simulation
# num_iter = 500
# mass_MC = []
# newa = []
# newt = []

# for b in range(num_iter):
    
#     A_err = np.random.normal(-3,3)
#     T_err = np.random.normal(-dT,dT)
#     new_Teff_p = teff_planet+T_err
#     new_a_p = age_planet+A_err
#     # print("L:\n", L_array)
#     # print("T:\n", T_array)
#     # print("New L: ", new_L_star, "New T:", new_Teff_star)
#     dist_calc = min_distance(new_a_p, new_Teff_p)
    
#     newa.append(new_a_p)
#     newt.append(new_Teff_p)  
#     mass_MC.append(dist_calc)
    
#     print("Iteration ", b+1, "out of ", num_iter)
    
# #%%
# #plots the corner plots
# samples = np.vstack([newa, newt, mass_MC])
# plt.rcParams['text.usetex'] = True
# plt.rcParams.update({'font.size': 22})
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'


# figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),labels = [r"$Age [Myr]$", r"$T_{eff} [K]$",r"$Mass [M_{J}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
# axes = np.array(figure.axes).reshape((3,3))


#%%

# def interp(df: pd.DataFrame):
#     temp = df['Teff(K)']
#     age = df['age(Gyr)']
#     mass = df['M/Msun']

#     x = np.array(age)
#     y = np.array(temp)
#     z = np.array(mass)
    
#     # isna = pd.isna(mass)
    
#     return interpolate.Rbf(x, y, z)

# #%%
# temp = data['Teff(K)']
# age = data['age(Gyr)']
# mass = data['M/Msun']

# x = np.array(age)
# y = np.array(temp)
# z = np.array(mass)
# mass_func = interpolate.Rbf(x,y,z)

# # mass_func(0.023,760)
# #%%
# # data = pd.read_csv("mass_evo_edit.txt", sep='\s+')
# f = interp(data)
# f(0.001,631)
#%%
# from numpy import savetxt
# MC_T = np.array(new_temps)
# MC_age = np.array(age_yr)
# MC_mass = np.array(mass_MC)

# dat_MC = np.array([MC_L, MC_T, MC_age, MC_mass])
# dat_MC = dat_MC.T
# savetxt('PARSEC_MC_0627.csv', dat_MC, delimiter = ',')

#%%
# mass_check = []

# for g in range(len(x)):
#     m_check = f(x[g],y[g])
#     mass_check.append(m_check)
# #%%
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'


# for q in range(0,15,1):
#     plt.plot(new_a[q], new_T[q], '-s', label = f"{mass_J[q]:.3}" r'$M_{J}$', markersize = 3, linewidth = 1.0)
#     plt.xlabel(r'$Age [Myr]$', fontsize = 30)
#     plt.ylabel(r'$Teff [K]$', fontsize = 30)
#     plt.axis([0, 50, 0, 1000])

# plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.5], [3.0]]), yerr = dT, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = '51 Eri b')
# plt.plot((x*1e9)/(1e6),y,'k.')
# plt.legend()
# plt.show()

# #%%
# #Plotting the base mass evolution tracks
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'

        
# for r in range(0,10,1):
#     plt.plot(age[r], T[r], '-s', label = f"{mass_J[r]:.2}" r'$M_{J}$', markersize = 3, linewidth = 1.0)
#     plt.xlabel(r'$Age [Myr]$', fontsize = 30)
#     plt.ylabel(r'$Teff [K]$', fontsize = 30)
#     # plt.axis([0, 50, 0, 1000])
# #plt.axis([3.9, 3.8,0.5,1.0])
# plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.5], [3.0]]), yerr = dT, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = '51 Eri b')
# plt.plot((x*1e9)/(1e6),y,'k.')
# plt.legend()
# plt.show()
    
# #%%
# #Plotting the base mass evolution tracks
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['xtick.direction'] = 'in'
# plt.rcParams['ytick.direction'] = 'in'

        
# # for r in range(5,5):
# plt.plot(age[5], T[5], '-s', label = f"{mass_J[5]:.2}" r'$M_{J}$', markersize = 3, linewidth = 1.0)
#     # plt.plot([r], T_c[r], '-s', label = f"{mass_J[r]:.2}" r'$M_{J}$ interp', markersize = 3, linewidth = 1.0)
# plt.xlabel(r'$Age [Myr]$', fontsize = 30)
# plt.ylabel(r'$Teff [K]$', fontsize = 30)
# plt.axis([0, 50, 0, 2000])
# #plt.axis([3.9, 3.8,0.5,1.0])
# plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.5], [3.0]]), yerr = dT, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = '51 Eri b')
# plt.plot((x*1e9)/(1e6),y,'k.')
# plt.legend()
# plt.show()
