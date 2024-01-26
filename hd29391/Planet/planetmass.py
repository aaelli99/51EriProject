import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator
from scipy.interpolate import interp2d
import pandas as pd
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Planet")
import numpy as np
import corner
#%%
#Code block to sort the table into individual mass tracks for plotting later
#Determines the beginning and ending indx position of each mass and sorts the 
#date using the indexes
data = pd.read_csv("mass_evo_edit.txt", sep='\s+')
Age = (data['age(Gyr)']*(1e9))/(1e6)
mass = []
age = []
logL = []
T = []
#masses from the data table
m = [0.0005,0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.006,
     0.007,0.008,0.009,0.01,0.011,0.0115,0.012,0.0125,0.013,0.0135,0.014,
     0.0145,0.015,0.0155,0.016,0.0165,0.017,0.0175,0.018,0.0185,0.019,0.0195,
     0.02,0.022,0.024,0.026,0.028,0.03,0.033,0.035,0.038,0.04,0.043,0.045,
     .048,0.05,0.053,0.055,0.058,0.06,0.063,0.064,0.065,0.066,0.067,0.068,
     0.069,0.07,0.071,0.072,0.073,0.074,0.075,0.077,0.078,0.08]
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
#Planet values
age_planet = 23.2                         #Age is in MYr
da_up = 1.72                             #Uncertainty is in Myr
da_down = 1.58
teff_planet = 807                       #Teff is in Kelvin
dT = 45
#loop converts the mass values into Jupiter masses instead of solar masses
mass_J = []
for i in range(len(m)):                 
    M = (m[i]*(1.989e33))/(1.899e30)
    mass_J.append(M)
    
    
#%%
#Testing out a new data point
age_planet = 1000
da_up = 15
da_down = 11
teff_planet = 900
dT = 48

#%%
#Code block to createa finer spacing of data points (not entirely useful 
#but helps visualize the curves better)
new_a = []
new_t = []
new_m = []
for k in range(0,len(m),1):
    mass_func = interp1d(age[k],T[k], kind = 'cubic')
    new_age = np.linspace(min(age[k]), max(age[k]), 10000)
    new_temps = mass_func(new_age)
    new_mass = np.ones(10000)*m[k]
    new_t.append(new_temps)
    new_a.append(new_age)
    new_m.append(new_mass)

#%%
new_a = []
new_t = []
new_m = []
for k in range(0,len(m),1):
    mass_func = interp1d(age[k],T[k], kind = 'linear')
    new_age = np.linspace(min(age[k]), max(age[k]), 10000)
    new_temps = mass_func(new_age)
    new_mass = np.ones(10000)*m[k]
    new_t.append(new_temps)
    new_a.append(new_age)
    new_m.append(new_mass)
#%%
#plotting the interpolated mass tracks
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.style.use('tableau-colorblind10')
# plt.figure(1)
plt.figure(figsize=(10, 7.5))

num_curves = 12
markers = ['.', 'v','o', 's', 'p', '8', 'P', 'D','<', 'X', 'h','>']
# cmap = plt.colormaps['cividis']
# colors = cmap(np.linspace(0,1,num_curves))
for q in range(2,num_curves,1):
    plt.plot(new_a[q], new_t[q],'-', label = f"{mass_J[q]:.2}" r'$\rm ~M_{Jup}$', 
             marker = markers[q],markersize = 7, linewidth = 1.0)
    # plt.plot(new_al[q], new_tl[q],'o--', label = f"{mass_J[q]:.3}" r'$M_{J} linear$', 
    #          markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$\rm Age~[Myr]$', fontsize = 30)
    plt.ylabel(r'$T_{\rm eff} \rm~[K]$', fontsize = 30)
    plt.axis([0, 50, 375, 1000])

plt.errorbar(age_planet, teff_planet, xerr = np.array([[1.58], [1.72]]), 
             yerr = dT, linestyle = 'None', linewidth = 2, capsize = 3, 
             color = 'k', label = r'51 Eri b')
# plt.grid(which = 'major')
# plt.grid(which='major', color='#DDDDDD', linewidth=0.8)
# plt.grid(which='minor', color='#EEEEEE', linestyle='-', linewidth=0.5)
plt.minorticks_on()
plt.legend()
plt.savefig('planetcurves_1220.pdf')
plt.show()
#%%
#Function to create an assymetrical random number generator based off of a 
#normal distribution
#Takes a random number generated from a mean of 0 and std of 1. If negative, 
#the lower uncertainty is used. 
#If positive the upper uncertainty is used
#The new uncertainty is then added to the actual value itself
#function returns the new_value
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
#using an interpolator to interpolate between the given curves

x = np.concatenate(new_a)
y = np.concatenate(new_t)
z = np.concatenate(new_m)


f = LinearNDInterpolator(list(zip(x,y)), z)

#%%
num_iter = 5000
mass_MC = []
newa = []
newt = []

for b in range(num_iter):
    A_err = assym_dist(da_down, da_up, age_planet)
    T_err = np.random.normal(teff_planet,dT)

    new_p_mass = f(A_err, T_err)
    print("Iteration ", b+1, "out of ", num_iter)
    
    newa.append(A_err)
    newt.append(T_err)  
    mass_MC.append(((new_p_mass)*(1.989e33))/(1.899e30))
    
#%%
#plotting the corner plot
samples = np.vstack([newa, newt, mass_MC])
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 25})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'


figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),
                       labels = [r"$\rm Age~[Myr]$", r"$T_{\rm eff}~[K]$",r"$\rm Mass$ $[M_{\rm Jup}]$"], 
                       show_titles = True, title_fmt = ".1f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((3,3))
ageq = corner.quantile(samples[0], q = (0.16, 0.5,0.84))
massq = corner.quantile(samples[2], q = (0.16, 0.5,0.84))
teffq = corner.quantile(samples[1], q = (0.16, 0.5,0.84))

ageup = ageq[2]-ageq[1]
agedown = ageq[1]-ageq[0]
massup = massq[2]-massq[1]
massdown = massq[1]-massq[0]
teffup = teffq[2]-teffq[1]
teffdown = teffq[1]-teffq[0]


axes[0, 0].set_title(r"$\rm Age~[Myr] =  %.1f ^{+%.1f}_{- %.1f}$" % (ageq[1], ageup, agedown), fontsize = 16, pad=10) 
axes[2, 2].set_title(r"$\rm Mass~[M_{\rm Jup}] =  %.1f ^{+%.1f}_{- %.1f}$" % (massq[1], massup, massdown),fontsize = 16, pad = 10)
axes[1, 1].set_title(r"$T_{\rm eff}~[K] =  %d ^{+%d}_{- %d}$" % (teffq[1], teffup, teffdown),fontsize = 16, pad = 10)


axes[2, 0].set_xticks([18,20,22,24,26,28])
axes[2, 1].set_xticks([650,700,750,800,850,900,950])
axes[2, 2].set_xticks([3.0,3.5,4.0,4.5,5.0,5.5,6.0])
axes[1, 0].set_xticks([18,20,22,24,26,28])
axes[1, 0].xaxis.set_tick_params(labelbottom=False)
axes[0, 0].set_xticks([18,20,22,24,26,28])
axes[0, 0].xaxis.set_tick_params(labelbottom=False)
axes[1, 1].set_xticks([650,700,750,800,850,900,950])
axes[1, 1].xaxis.set_tick_params(labelbottom=False)

axes[1, 0].set_yticks([650,700,750,800,850,900,950])
axes[2, 0].set_yticks([3.0,3.5,4.0,4.5,5.0,5.5,6.0])
axes[2, 1].set_yticks([3.0,3.5,4.0,4.5,5.0,5.5,6.0])
axes[2, 1].yaxis.set_tick_params(labelleft=False)
# figure.title('LinearND Interpolation')
# figure.savefig('planet_corner_1005.pdf')
#%%
#function that does a one dimensional interpolation of the age and temperature 
#arrays from the mass evolution table
#It interpolates then a new temperature is determined off of the function the 
#interpolation created when the target age
#is inputted. 
#outputs the mass array and the temperature array
def interp_t(age,T, target_Age):
    interpt = []
    masses = []
    for i in range(len(age)):
        interp_t = interp1d(age[i],T[i], kind = 'linear')

        if target_Age>age[i][0]:       
            new_T = interp_t(target_Age)
            masses.append(m[i])
            interpt.append(new_T)
    M = np.array(masses)
    t = np.array(interpt)
    
    return((t,M))

#%%
#Monte Carlo simulation to obtain uncertainties
#A_err is the new age value based off of the asymmetrical normal distribution 
#function above
#T_err is the new teff value based off of the random.normal distribution with 
#teff_planet as the mean and dT as the sigma
#used the interpolation function above to get the new temperature and mass arrays
#then performs a one dimensional interpolation between the new temperature 
#array and the mass array created from the 
#target age
#then the random teff value is plugged into the new function and read out to 
#get a new mass.
#values are appended to lists for plotting in a corner plot later on
niter = 1000
mass_MC = []
newa = []
newt = []

for i in range(niter):
    A_err = assym_dist(da_down, da_up, age_planet)
    T_err = np.random.normal(teff_planet,dT)
    # print("Age:", A_err)
    # print("Temp:", T_err)
    temp_func = interp_t(age,T,A_err)
    temps = temp_func[0]
    mm = temp_func[1]
    mass_func = interp1d(temps, mm, kind = 'linear')
    newm = mass_func(T_err)
    # print('Planet mass:', (newm*(1.989e30))/(1.899e27))
    print("Iteration ", i+1, "out of ", niter)
    newa.append(A_err)
    newt.append(T_err)
    mass_MC.append((newm*(1.989e30))/(1.899e27))


        
        
        