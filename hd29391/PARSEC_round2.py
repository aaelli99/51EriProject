#Re-do of the PARSEC isochrone codes
#Ashley Elliott
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate
from scipy.interpolate import interp1d
import random
import corner
#%%
#Removing the unnecessary lines and combining to fit one file
#ONLY run if using a new set of data
#f = open("PARSEC_data_66_88_new.txt", 'r')
f = open('PARSEC_data_0426.txt','r')
#f_new = open("PARSEC_data_edit_66_88_new.txt", 'a')
f_new = open("PARSEC_data_edit_0426.txt",'a')
f_new.write("Zini MH logAge Mini int_IMF Mass logL logTe logg label McoreTP C_O period0  period1  period2  period3  period4  pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  Umag    Bmag    Vmag    Rmag    Imag    Jmag    Hmag    Kmag \n")
no = "#"
for line in f:
    if not no in line:
        f_new.write(line)
        f_new.write("\n")

f_new.close()      
f.close()
#%%
#Cell reads in the data and determines the beginning and ending indices for each age in order to sort the data into separate arrays based on the age it goes with
#data=pd.read_csv("PARSEC_data_edit_66_88_new.txt", sep='\s+')
data = pd.read_csv("PARSEC_data_edit_0426.txt", sep='\s+')
#removing columns I wont need
data.drop(['label','McoreTP','C_O','period0','period1', 'period2','period3','period4','pmode','Mloss', 'tau1m','X','Y','Xc','Xn','Xo','Cexcess','Z','mbolmag','Umag', 'Bmag','Vmag','Rmag','Imag','Jmag','Hmag'],axis=1, inplace=True)
logL = []
logTeff = []
mass = []
#age = [6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6]
#age = np.round(np.log10(np.arange(1e7,5.2e7,2e6)),5)
age = np.around(np.arange(6,8,0.05),3)
lens = np.empty(len(age))

for ii in range(len(lens)):
    idx = 0
    #print("Age Array:", age[ii])
    for i in range(len(data['logAge'])):
        if round(data["logAge"][i],2) == age[ii]:
            idx = idx+1
        lens[ii] = idx
    beg_pos = np.empty(len(age))
    end_pos = np.empty(len(age))
    tot = lens[0]-1
    beg_pos[0] = 0
    end_pos[0] = int(tot)
    for x in range(len(lens)-1):
        beg_pos[x+1] = int(tot+1)
        tot = tot + lens[x+1]
        end_pos[x+1] = int(tot)
    
    # print("Begin: ",beg_pos[ii])
    # print("End: ", end_pos[ii]) 

#%%
for y in range(len(age)):
    lums = data["logL"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    temps = data["logTe"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    Mass = data["Mass"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    logL.append(lums)
    logTeff.append(temps)
    mass.append(Mass)

#%%
#Star information
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)

#%%
#Applying luminosity cuts to the data more than 2.5 luminosity
new_L_cut = []
new_Teff_cut = []
new_mass_cut = []

for i in range(0,len(age),1):
    lum_full = logL[i]
    teff_full = logTeff[i]
    mass_full = mass[i]
    pairs = list(zip(lum_full, teff_full, mass_full))
    new_pairs = []
    
    for i in range(len(lum_full)):
        if pairs[i][0] <= 2.5 and pairs[i][0]>-9.999:
            # print("It worked")
            # print("Lum is:", pairs[i][0])
            new_pairs.append((pairs[i][0], pairs[i][1], pairs[i][2]))

    logL_cut, logTeff_cut, mass_cut = zip(*new_pairs)
    
    new_L_cut.append(logL_cut)
    new_Teff_cut.append(logTeff_cut)
    new_mass_cut.append(mass_cut)
    
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for jj in range(len(age)):
    plt.plot(new_Teff_cut[jj], new_L_cut[jj], '-', label = f"{10**(age[jj])}", markersize = 5.0, linewidth = 0.2)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    # plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = 'HD29391')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()
#%%

new_M = []
new_T = []
new_L = []

for z in range(len(age)):
    newt = np.linspace(min(new_Teff_cut[z]),max(new_Teff_cut[z]), 10000)
    f_Age = interp1d(new_Teff_cut[z], new_L_cut[z], kind = 'linear')
    f_Mass = interp1d(new_Teff_cut[z], new_mass_cut[z], kind = 'linear')
    
    newl = np.empty(len(newt))
    newm = np.empty(len(newt))
    lum_int = f_Age(newt)
    mass_int = f_Mass(newt)
    
    new_M.append(mass_int)
    new_L.append(lum_int)
    new_T.append(newt)
    
#%%

new_M = []
new_T = []
new_L = []

for z in range(len(age)):
    newm = np.linspace(min(new_mass_cut[z]),max(new_mass_cut[z]), 10000)
    f_Lum = interp1d(new_mass_cut[z], new_L_cut[z], kind = 'linear')
    f_Temp = interp1d(new_mass_cut[z], new_Teff_cut[z], kind = 'linear')
    
    newl = np.empty(len(newm))
    newt = np.empty(len(newm))
    lum_int = f_Lum(newm)
    temp_int = f_Temp(newm)
    
    new_M.append(newm)
    new_L.append(lum_int)
    new_T.append(temp_int)

#%%
#minimum distance function

def min_distance(x1,y1,T_array,L_array):
    dist = []
    min_dist = []     
    for m in range(len(age)):
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
    
    print("The minimum distance between points is", min_iso)
    print("Age of isochrone:",f"{10**(age[pos_iso]):.2e} years")
    return((age_iso, mass_iso_L, mass_iso_T, L, T, M))

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for jj in range(len(age)):
    plt.plot(new_T[jj], new_L[jj], '.', label = f"{round((10**(age[jj]))/10**6,1)} Myr", markersize = 5.0, linewidth = 0.2)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = 'HD29391')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()
#%%
niter = 1000
age_MC = []
massL_MC = []
massT_MC = []
newvals = []

for ii in range(niter):
    L_err = np.random.normal(-dL,dL)
    T_err = np.random.normal(-dT,dT)
    
    new_Teff_star = Teff_star+T_err
    new_L_star = L_star+L_err
    
    dist_calc = min_distance(new_Teff_star, new_L_star, new_T, new_L)
    
    newvals.append((new_Teff_star,new_L_star))
    
    age_MC.append(dist_calc[0])
    massL_MC.append(dist_calc[1])
    massT_MC.append(dist_calc[2])
    
    print("Iteration ", ii, "out of ", niter)
    
#%%
new_temps = []
age_yr = []
lum_nl = []
for j in range(len(age_MC)):
    t = 10**(newvals[j][0])
    a = (10**(age_MC[j]))/(1e6)
    l = 10**(newvals[j][1])
    new_temps.append(t)
    age_yr.append(a)
    lum_nl.append(l)
    
    
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

samples = np.vstack([age_yr, massL_MC, new_temps, lum_nl])

figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),range = [(18,30),(1.53,1.55), (7320,7350),(5.6,5.8)],labels = [r"$Age [Myr]$", r"$Mass [M_{\odot}]$", r"$T_{eff} [K]$", r"$L_{\star} [L_{\odot}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((4,4))

value1 = np.mean(samples.T,axis=0)
    
    