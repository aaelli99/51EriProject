#Ashley Elliott
#PARSEC isochrone data for HD 29391
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
#import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
import random
import corner
from scipy.spatial import cKDTree

# from isochrones.interp import DFInterpolator as DF
#%%
#Removing the unnecessary lines and combining to fit one file
#ONLY run if using a new set of data
#f = open("PARSEC_data_66_88_new.txt", 'r')
#f = open('PARSEC_data_1_5.txt','r')
f = open('PARSEC_data_0620.txt','r')
# f = open('PARSEC_data_0614_6585.txt')
#f_new = open("PARSEC_data_edit_66_88_new.txt", 'a')
#f_new = open("PARSEC_data_edit_1_5.txt",'a')
# f_new = open('PARSEC_data_0614_6585_edit.txt','a')
f_new = open('PARSEC_data_0620_edit.txt','a')
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
# data = pd.read_csv("PARSEC_data_0614_6585_edit.txt", sep='\s+')
data = pd.read_csv("PARSEC_data_0620_edit.txt", sep='\s+')
# data1 = pd.read_csv("PARSEC_data_0616_edit.txt", sep='\s+')
#removing columns I wont need
data.drop(['label','McoreTP','C_O','period0','period1', 'period2','period3','period4','pmode','Mloss', 'tau1m','X','Y','Xc','Xn','Xo','Cexcess','Z','mbolmag','Umag', 'Bmag','Vmag','Rmag','Imag','Jmag','Hmag'],axis=1, inplace=True)
# data1.drop(['label','McoreTP','C_O','period0','period1', 'period2','period3','period4','pmode','Mloss', 'tau1m','X','Y','Xc','Xn','Xo','Cexcess','Z','mbolmag','Umag', 'Bmag','Vmag','Rmag','Imag','Jmag','Hmag'],axis=1, inplace=True)
logL = []
logTeff = []
mass = []
# logL1 = []
# logTeff1 = []
# mass1 = []
#age = [6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6]
# age = np.round(np.log10(np.arange(1e7,5.2e7,2e6)),5)
age = np.round(np.log10(np.arange(1e7,5.1e7,1e6)),5)
# age = np.arange(6.5,7.65,0.05)
# age1 = np.arange(6.5,8.6,0.1)
#age = 6.6
lens = np.empty(len(age))
# lens1 = np.empty(len(age1))
for ii in range(len(age)):
    idx = 0
    for i in range(len(data['logL'])):
        if data["logAge"][i] == age[ii]:
            # print("Data is: ", data["logAge"][i], "with index ", i)
            # print("Age is: ", age[ii])
            idx = idx+1
            #print(idx)
        lens[ii] = idx
    
    beg_pos = np.empty(len(age))
    end_pos = np.empty(len(age))
    # print(lens[ii])
    tot = lens[0]-1
    beg_pos[0] = 0
    end_pos[0] = int(tot)
    for x in range(len(lens)-1):
        # print(x)
        beg_pos[x+1] = int(tot+1)
        tot = tot + lens[x+1]
        end_pos[x+1] = int(tot)
        # print(tot)   
        
    # print("Begin: ",beg_pos[ii])
    # print("End: ", end_pos[ii]) 
# for jj in range(len(age1)):
#     idx1 = 0
#     for j in range(len(data1['logL'])):
#         if np.round(data1["logAge"][j],2) == np.round(age1[jj],2):
#             # print("Data is: ", data["logAge"][i], "with index ", i)
#             # print("Age is: ", age[ii])
#             idx1 = idx1+1
#             #print(idx)
#         lens1[jj] = idx1
    
#     beg_pos1 = np.empty(len(age1))
#     end_pos1 = np.empty(len(age1))
#     #print(lens[ii])
#     tot1 = lens1[0]-1
#     beg_pos1[0] = 0
#     end_pos1[0] = int(tot1)
#     for w in range(len(lens1)-1):
#         # print(x)
#         beg_pos1[w+1] = int(tot1+1)
#         tot1 = tot1 + lens1[w+1]
#         end_pos1[w+1] = int(tot1)
#         # print(tot)   
        
#     # print("Begin: ",beg_pos[ii])
#     # print("End: ", end_pos[ii]) 
        
    
for y in range(len(age)):
#     print("Index:", y)
#     print("Begin: ",beg_pos[y])
#     print("End: ", end_pos[y])

    lums = data["logL"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    temps = data["logTe"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    Mass = data["Mass"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    logL.append(lums)
    logTeff.append(temps)
    mass.append(Mass)
    
# for z in range(len(age1)):
# #     print("Index:", y) 
# #     print("Begin: ",beg_pos[y])
# #     print("End: ", end_pos[y])

#     lums1 = data1["logL"][int(beg_pos1[z]):int(end_pos1[z])+1].values.tolist()
#     temps1 = data1["logTe"][int(beg_pos1[z]):int(end_pos1[z])+1].values.tolist()
#     Mass1 = data1["Mass"][int(beg_pos1[z]):int(end_pos1[z])+1].values.tolist()
#     logL1.append(lums1)
#     logTeff1.append(temps1)
#     mass1.append(Mass1)
    
    
#%%
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)

for x in range(0,24,1):
    plt.plot(logTeff[x], logL[x], '-s', label = f"{10**(age[x]):.3e} years", markersize = 3, linewidth = 1.0)
    # plt.plot(logTeff1[x], logL1[x],'-s', label = f"{10**(age1[x]):.3e} years (new data)", markersize = 3, linewidth = 1.0)
    # plt.plot(logTeff[5], logL[5], '-g', label = f"{10**(age[5]):.3e} years (1_5 data)", markersize = 3, linewidth = 1.0)
    # plt.plot(logTeff1[8], logL1[8],'-c', label = f"{10**(age1[8]):.3e} years (new data)", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([4.5, 3.3, -2.5, 6])
  #plt.axis([3.9, 3.8,0.5,1.0])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()
#%%

#df = pd.DataFrame()

#df[f"{10**(age[5]):.2e} years"] = logTeff[5]

#df.plot()
I = []
for i in range(0,20,1):
    idx = np.empty(len(logTeff[i]))
    for ii in range(len(logTeff[i])):
        idx[ii] = ii
    I.append(idx)
    
for k in range(0,20,1):
    plt.plot(I[k], logTeff[k],'-', label = f"{10**(age[k]):.2e} years")   
    plt.xlabel('Index')
    plt.ylabel(r'$log(T_{eff})$', fontsize = 30)
plt.legend()
plt.show()
        
#%%
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)
#%%
#plots the isochrones based on data
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
L_star = np.log10(5.712)
Teff_star = np.log10(7414)
dL = 0.434*(0.096/5.712)
dT = 0.434*(32/7414)

for x in range(0,20,1):
    plt.plot(logTeff[x], logL[x], '-s', label = f"{10**(age[x]):.2e} years", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    #plt.axis([3.88, 3.86, 0.6, 0.9])
    #plt.axis([3.9, 3.8,0.5,1.0])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%
plt.plot(logTeff[9], logL[9], '-s', label = f"{10**(age[9]):.2e} years", markersize = 5, linewidth = 1.0)
plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
plt.axis([3.88, 3.86, 0.6, 0.9])
#plt.axis([3.9, 3.8,0.5,1.0])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%
#calculates the distance between my data point and the points in the curves
x1 = L_star
y1 = Teff_star
dist = []
min_dist = []
for m in range(0,21,1):
    chrone_L = logL[m]
    chrone_T = logTeff[m]
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
#written 02/01/2022
#cutting the luminosity values to any less than 2.5 because that's all we care about for now due to the isochrones entering post-main sequence

new_L_cut = []
new_Teff_cut = []
new_mass_cut = []

for i in range(len(age)):
    lum_full = logL[i]
    teff_full = logTeff[i]
    mass_full = mass[i]
    pairs = list(zip(lum_full, teff_full, mass_full))
    new_pairs = []
    
    for i in range(len(lum_full)):
        if pairs[i][0] <= 2.5 and pairs[i][0]>-9.999:
        #print("It worked")
        #print("Lum is:", pairs[i][0])
            new_pairs.append((pairs[i][0], pairs[i][1], pairs[i][2]))

    logL_cut, logTeff_cut, mass_cut = zip(*new_pairs)
    
    new_L_cut.append(logL_cut)
    new_Teff_cut.append(logTeff_cut)
    new_mass_cut.append(mass_cut)

    

#%%
#interpolates the data so we get a finer grid of data points
new_L = []
new_T = []
new_M = []


for z in range(0,len(age),1):
    age_func = interp1d(new_Teff_cut[z], new_L_cut[z], kind = 'linear')
    mass_func = interp1d(new_Teff_cut[z], new_mass_cut[z])
    newtemps = np.arange(min(new_Teff_cut[z]), max(new_Teff_cut[z])-0.0001, 0.0001)
    #newtemps = np.linspace(min(new_Teff_cut[z]), max(new_Teff_cut[z]), 100)
    newlum = np.empty(len(newtemps))
    new_mass = np.empty(len(newtemps))
    new_mass = mass_func(newtemps)
    newlum = age_func(newtemps)
    new_L.append(newlum)
    new_T.append(newtemps)
    new_M.append(new_mass)
   # print("LogTeff:", logTeff[g])
#%%
new_L = []
new_T = []
new_M = []

for g in range(1):  
    #interp = RegularGridInterpolator((age[g], mass[g]), (logL[g], logTeff[g]) )
    age_func = interp1d(logTeff[g], logL[g], kind = 'linear')
    mass_func = interp1d(logTeff[g], mass[g])
    #newtemps = np.arange(min(logTeff[g]), max(logTeff[g]), 0.0001)
    newtemps = np.linspace(max(logTeff[g]), min(logTeff[g]), 100)
    newlum = np.empty(len(newtemps))
    new_mass = np.empty(len(newtemps))
    new_mass = mass_func(newtemps)
    newlum = age_func(newtemps)
    new_L.append(newlum)
    new_T.append(newtemps)
    new_M.append(new_mass)
   # print("LogTeff:", logTeff[g])
    #print(logL[0])

    
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for jj in range(0,15,1):
    plt.plot(new_T[jj], new_L[jj], '.', label = f"{round((10**(age[jj]))/10**6,1)} Myr", markersize = 5.0, linewidth = 0.2)
    #plt.plot(logTeff[jj], logL[jj], '-s', label = f"{10**(age[jj]):.2e} years (Noninterp)", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = 'HD29391')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()


#%%
#Rewriting the minimum distance function
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
        
    print("The closest isochrone is:", 10**(min_age), "Myr") 
    # print("The minimum distance is:", min_iso) 
    # print("The index is:", min_pos)
    print("The mass is:", min_mass, "solar masses")
    
    return((min_age,min_mass))

#%%
    # #snippet to find the closest mass
    
    # newm = np.linspace(min(new_mass_cut[min_pos]),max(new_mass_cut[min_pos]), 10000)
    # f_Lum = interp1d(new_mass_cut[min_pos], new_L_cut[min_pos], kind = 'linear')
    # f_Temp = interp1d(new_mass_cut[min_pos], new_Teff_cut[min_pos], kind = 'linear')
        
    # # newl = np.empty(len(newm))
    # # newt = np.empty(len(newm))
    # lum_int = f_Lum(newm)
    # temp_int = f_Temp(newm)
    
    # mass_M = newm
    # mass_L = lum_int
    # mass_T = temp_int
    
    # # new_vals_m = list(zip(mass_M,mass_L,mass_T))
    
    # #ChatGPT code
    # points = np.column_stack((mass_T, mass_L))
    # tree = cKDTree(points)
    
    # # Find the indices of the nearest neighbors to the specified teff and luminosity values
    # dist, indices = tree.query((Teff_star, L_star), k=1)
    
    # # Retrieve the corresponding age and mass values
    # interp_mass = mass_M[indices]
    
    # print("The mass of the star is:", interp_mass)
    
    # return((min_age,interp_mass))
    

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

plt.plot(mass_T[0], mass_L[0],'.', label = f"{round((10**(age[min_pos]))/10**6,1)} Myr", markersize = 5.0, linewidth = 0.2)
    # plt.plot(logTeff[jj], logL[jj], '-s', label = f"{10**(age[jj]):.2e} years (Noninterp)", markersize = 3, linewidth = 1.0)
plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
#plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = 'HD29391')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()

#%%
#minimum distance function

# def min_distance(x1,y1,T_array,L_array):
#     dist = []
#     min_dist = []
#     # chrone_L = L_array
#     # chrone_T = T_array
#     # d = np.empty(len(chrone_L))
#     # for k in range(len(chrone_L)):
#     #     x2 = chrone_T[k]
#     #     y2 = chrone_L[k]
#     #     d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
#     # dist.append(d)
                     
#     for m in range(0,len(age),1):
#         chrone_L = L_array[m]
#         chrone_T = T_array[m]
#         d = np.empty(len(chrone_L))
#         for k in range(len(chrone_L)):
#             y2 = chrone_L[k]
#             x2 = chrone_T[k]
#             d[k] = np.sqrt((abs(x1-x2)**2)+(abs(y1-y2))**2)
    
#         dist.append(d)
    
#     for n in range(len(dist)):
#         min_dist.append(np.min(dist[n]))

#     min_iso = np.min(min_dist)
#     pos_iso = np.argmin(min_dist)

#     L = new_L[pos_iso]
#     T = new_T[pos_iso]
#     M = new_M[pos_iso]
    
#     l = np.asarray(L)
#     t = np.asarray(T)
#     idx_l = (np.abs(l-L_star)).argmin()
#     idx_t = (np.abs(t-Teff_star)).argmin()
    
#     age_iso = age[pos_iso]
#     mass_iso_L = M[idx_l]
#     mass_iso_T = M[idx_t]
    
#     return((age_iso, mass_iso_L, mass_iso_T, L, T, M))

#%%
#Monte Carlo simulations
num_iter = 1000
age_MC = []
mass_MC = []
# massT_MC = []
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
    # massT_MC.append(dist_calc[2])
    
    print("Iteration ", b+1, "out of ", num_iter)

#%%
new_temps = []
age_yr = []
lum_nl = []
for i in range(num_iter):
    t = 10**(newvals[i][0])
    a = (10**(age_MC[i]))/(1000000)
    l = 10**(newvals[i][1])
    new_temps.append(t)
    age_yr.append(a)
    lum_nl.append(l)
#%%
samples = np.vstack([age_yr, mass_MC, new_temps, lum_nl])

# mean_age = np.mean(age_MC)
# mean_mass = np.mean(massT_MC)
# std_age = np.std(age_MC)
# std_mass = np.std(massT_MC)

#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),range = [(10,30),(1.5,1.575), (7250,7500),(5.3,5.9)],labels = [r"$Age [Myr]$", r"$Mass [M_{\odot}]$", r"$T_{eff} [K]$", r"$L_{\star} [L_{\odot}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((4,4))

value1 = np.mean(samples.T,axis=0)

#%%
from numpy import savetxt
MC_L = np.array(lum_nl)
MC_T = np.array(new_temps)
MC_age = np.array(age_yr)
MC_mass = np.array(mass_MC)

dat_MC = np.array([MC_L, MC_T, MC_age, MC_mass])
dat_MC = dat_MC.T
savetxt('PARSEC_MC_0619_v4.csv', dat_MC, delimiter = ',')

#%%
y1 = L_star
x1 = Teff_star
dist = []
min_dist = []
for m in range(0,len(age),1):
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
#ChatGPT basic calculation
def star_age(temperature, luminosity):
    # Calculate the star's main sequence lifetime
    log_t_ms = 10.13 + (1.62 * np.log10(temperature/5770))
    t_ms = 10**log_t_ms
    # Calculate the star's current age
    age = t_ms * (1 - (luminosity/3.846e26)**(1/2.5))
    return age

# Example usage
temperature = 7414
luminosity = 5.712
print("Star age:", star_age(temperature, luminosity), "years")

#%%
#ChatGPT more advanced model
from isochrones import StarModel
def star_age(temperature, luminosity):
    # Load the MIST isochrones
    mist = StarModel('mist')
    # Convert temperature and luminosity to logarithmic units
    log_temperature = np.log10(temperature)
    log_luminosity = np.log10(luminosity)
    # Interpolate the isochrones to get an estimate of the star's age
    result = mist.interp_value([log_temperature, log_luminosity], ['log10_isochrone_age'])
    log_age = result['log10_isochrone_age'][0]
    age = 10**log_age
    return age
# Example usage
temperature = 7414
luminosity = 5.712
print("Star age:", star_age(temperature, luminosity), "years")

#%%
#Chatgpt more advance with PARSEC model
def star_age(temperature, luminosity):
    # Load the PARSEC isochrones
    parsec = StarModel('parsec')
    # Convert temperature and luminosity to logarithmic units
    log_temperature = np.log10(temperature)
    log_luminosity = np.log10(luminosity)
    # Interpolate the isochrones to get an estimate of the star's age
    result = parsec.interp_value([log_temperature, log_luminosity], ['log10_isochrone_age'])
    log_age = result['log10_isochrone_age'][0]
    age = 10**log_age
    return age
# Example usage
temperature = 7414
luminosity = 5.712
print("Star age:", star_age(temperature, luminosity), "years")

#%%
#ChatGPT interpolation for mass
def star_properties(temperature, luminosity):
    # Load the PARSEC isochrones
    parsec = StarModel('parsec')
    # Convert temperature and luminosity to logarithmic units
    log_temperature = np.log10(temperature)
    log_luminosity = np.log10(luminosity)
    # Interpolate the isochrones to get estimates of the star's age and mass
    result = parsec.interp_value([log_temperature, log_luminosity], ['log10_isochrone_age', 'mass'])
    log_age = result['log10_isochrone_age'][0]
    age = 10**log_age
    mass = result['mass'][0]
    return age, mass
# Example usage
temperature = 5778
luminosity = 3.846e26
age, mass = star_properties(temperature, luminosity)
print("Star age:", age, "years")
print("Star mass:", mass, "solar masses")