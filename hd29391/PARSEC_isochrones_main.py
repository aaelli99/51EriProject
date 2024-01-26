#PARSEC Isochrone code
#############################################################################
#     Code to read in PARSEC isochrones, interpolate, and to run Monte      #
#     Carlo simulations to get uncertainties.                               #
#     Steps:                                                                #
#     1. Read in PARSEC isochrone file                                      #
#     2. Edit the file to get the data I want (removes unnecessary columns) #
#     3. Extract each year's isochrone by finding the beginning and end     #
#        index and then sorting the data by age                             #
#     4. Reads in the star information                                      #
#     5. Plots to check the isochrones read in properly and isochrones      #
#        look correct                                                       #
#     6. Cuts the data to prevent any luminosity over 2.5                   #
#     7. Interpolates the isochrone lines to create a smaller step size     #
#        for each line (i.e. 10 points in a line to 1000 points, etc.)      #
#     8. Plots new interpolated isochrones                                  #
#     9. Distance minimization function that reads in the interpolated      #
#        isochrones and calculates the distance between the star point and  #
#        each line on the curve. Finds the minimum distance and stores it.  #
#        Then finds the isochrone with the smallest distance, stores the    #
#        point and its position. Finds the corresponding mass value for     #
#        that point. Returns the age and the mass                           #
#    10. Monte Carlo run that varies the star's luminosity and temp by      #
#        some random value in the uncertainty range and runs the distance   #
#        minimization function on it. Stores the age, mass, temp, and lum.  #
#    11. Creates a corner plot that shows the distribution and returns the  #
#        uncertainties for each value.                                      #
#                                                                           #
#############################################################################


#%%
#Importing all necessary packages
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\Isochrone Data")
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy import interpolate
import pandas as pd
from scipy.interpolate import interp1d
import random
import corner
from scipy.interpolate import LinearNDInterpolator
import matplotlib.ticker as ticker
#%%
#Reading in new data file and rewriting to new file that contains only the data
#ONLY to be used once per new data file otherwise it will cause issues
f = open('PARSEC_data_0623_newZ.txt','r')            #Current isochrone file
#f_new = open("PARSEC_data_edit_1_5.txt",'a')  #Isochrone that used to work
f_new = open('PARSEC_data_0623_edit.txt','a')
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
data = pd.read_csv("PARSEC_data_0623_edit.txt", sep='\s+')
#removing columns I wont need
data.drop(['McoreTP','C_O','period0','period1', 'period2','period3','period4','pmode','Mloss', 'tau1m','X','Y','Xc','Xn','Xo','Cexcess','Z','mbolmag','Umag', 'Bmag','Vmag','Rmag','Imag','Jmag','Hmag'],axis=1, inplace=True)

logL = []
logTeff = []
mass = []
label = []
age = np.round(np.log10(np.arange(1e7,5.05e7,5e5)),5)
# age = np.arange(6.5,7.65,0.05)
lens = np.empty(len(age))
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
    
for y in range(len(age)):
#     print("Index:", y)
#     print("Begin: ",beg_pos[y])
#     print("End: ", end_pos[y])

    lums = data["logL"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    temps = data["logTe"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    Mass = data["Mass"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    lab = data["label"][int(beg_pos[y]):int(end_pos[y])+1].values.tolist()
    logL.append(lums)
    logTeff.append(temps)
    mass.append(Mass)
    label.append(lab)

#%%
#plotting the isochrones before interpolation
plt.rcParams.update({'font.size': 15})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
L_star = np.log10(5.712)
Teff_star = np.log10(7424)
dL = 0.434*(0.096/5.712)
dT = 0.434*(45/7424)

for x in range(0,len(age),1):
    plt.plot(logTeff[x], logL[x], '-s', label = f"{10**(age[x]):.3e} years", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([4.5, 3.3, -2.5, 6])
plt.axis([3.9, 3.8,0.5,1.0])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k', label = 'HD29391')
plt.legend()
plt.show()

#%%
#Star information
#new star info date: 06/27
L_star = np.log10(5.723)
Teff_star = np.log10(7424)
dL = 0.434*(0.096/5.723)
dT = 0.434*(45/7424)
R_star = 1.45
dR = 0.01
t_star = 7424
dt_star = 45

#%%
#written 02/01/2022
#cutting the luminosity values to any less than 2.5 because that's all we care about for now due to the isochrones entering post-main sequence
#Added to the if statement for the larger than -9.999 to prevent issues

new_L_cut = []
new_Teff_cut = []
new_mass_cut = []
new_label_cut = []
for i in range(len(age)):
    lum_full = logL[i]
    teff_full = logTeff[i]
    mass_full = mass[i]
    lab_full = label[i]
    pairs = list(zip(lum_full, teff_full, mass_full, lab_full))
    new_pairs = []
    
    for i in range(len(lum_full)):
        if pairs[i][0] <= 2.5 and pairs[i][0]>-9.999:
        #print("It worked")
        #print("Lum is:", pairs[i][0])
            new_pairs.append((pairs[i][0], pairs[i][1], pairs[i][2], pairs[i][3]))

    logL_cut, logTeff_cut, mass_cut, lab_cut = zip(*new_pairs)
    
    new_L_cut.append(logL_cut)
    new_Teff_cut.append(logTeff_cut)
    new_mass_cut.append(mass_cut)
    new_label_cut.append(lab_cut)

#%%
sigma = 5.67e-8     #W/m^2K^4
rteff = []
rlum = []
rmass = []
radius = []
for dum in range(len(new_L_cut)):
    rT = []
    rL = []
    R = []
    masss = []
    for dumm in range(len(new_L_cut[dum])):
        r_teff = 10**(new_Teff_cut[dum][dumm])
        r_lum = 10**(new_L_cut[dum][dumm])*(3.827e26)
        rads = (np.sqrt(r_lum/(4*np.pi*sigma*(r_teff**4))))/(6.96e8)
        mas = new_mass_cut[dum][dumm]*(1.989e33)
        rT.append(r_teff)
        rL.append(r_lum)
        R.append(rads)
        masss.append(mas)
    rteff.append(rT)
    rlum.append(rL)
    rmass.append(masss)
    radius.append(R)
    

#%%
#Interpolates the isochrones to get a finer grid spacing for the lines
new_R = []
new_T = []
new_M = []
new_A = []

for z in range(0,len(age),1):
    age_func = interp1d(rteff[z], radius[z], kind = 'linear')
    mass_func = interp1d(radius[z], rmass[z], kind = 'linear')
    newtemps = np.arange(min(rteff[z]), max(rteff[z]), 1)
    newrads = np.empty(len(newtemps))
    new_a = np.ones(len(newtemps))*age[z]
    newrad = age_func(newtemps)
    new_mass = mass_func(newrad)
    new_R.append(newrad)
    new_T.append(newtemps)
    new_M.append(new_mass)
    new_A.append(new_a)
#%%
#Interpolates the isochrones to get a finer grid spacing for the lines
new_L = []
new_T = []
new_M = []
new_A = []

for z in range(0,len(age),1):
    age_func = interp1d(new_Teff_cut[z], new_L_cut[z], kind = 'linear')
    mass_func = interp1d(new_L_cut[z], new_mass_cut[z])
    newtemps = np.arange(min(new_Teff_cut[z]), max(new_Teff_cut[z])-0.0001, 0.0001)
    newlum = np.empty(len(newtemps))
    new_a = np.ones(len(newtemps))*age[z]
    newlum = age_func(newtemps)
    new_mass = mass_func(newlum)
    new_L.append(newlum)
    new_T.append(newtemps)
    new_M.append(new_mass)
    new_A.append(new_a)
#%%
# for z in range(0,len(age),1):
#     mass_func = interp1d(new_A[z], new_mass_cut[z])
#     newage = np.arange(min(new_A[z]), max(new_A[z]), 0.0001)
#     new_mass = mass_func(newage)
#     new_M.append(new_mass)

#%%
#plots the interpolated isochrones
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for jj in range(20,32,1):
    plt.plot(new_T[jj], new_R[jj], '-s', label = f"{round((10**(age[jj]))/10**6,1)} Myr", markersize = 5.0, linewidth = 0.2)
    #plt.plot(logTeff[jj], logL[jj], '-s', label = f"{10**(age[jj]):.2e} years (Noninterp)", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$T_{eff}$', fontsize = 30)
    plt.ylabel(r'$R/R_{\odot})$', fontsize = 30)
    plt.axis([8000, 6000, 1.00, 2.0])
plt.errorbar(t_star, R_star, yerr = dR, xerr = dt_star,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = '51 Eridani')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()
#%%
#plots the interpolated isochrones
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

for jj in range(20,32,1):
    plt.plot(new_T[jj], new_L[jj], '-s', label = f"{round((10**(age[jj]))/10**6,1)} Myr", markersize = 5.0, linewidth = 0.2)
    #plt.plot(logTeff[jj], logL[jj], '-s', label = f"{10**(age[jj]):.2e} years (Noninterp)", markersize = 3, linewidth = 1.0)
    plt.xlabel(r'$log(T_{eff})$', fontsize = 30)
    plt.ylabel(r'$log(L/L_{\odot})$', fontsize = 30)
    plt.axis([3.88, 3.86, 0.6, 0.9])
plt.errorbar(Teff_star, L_star, yerr = dL, xerr = dT,linestyle = 'None', linewidth = 3, capsize = 3, color = 'k', label = '51 Eridani')
plt.title('PARSEC Isochrones')
plt.legend(prop={'size': 15})
plt.show()


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
        
    print("The closest isochrone is:", 10**(min_age), "Myr") 
    # print("The minimum distance is:", min_iso) 
    # print("The index is:", min_pos)
    print("The mass is:", min_mass, "solar masses")
    
    return((min_age,min_mass))

#%%
x = np.concatenate(new_T)
y = np.concatenate(new_L)
z = np.concatenate(new_A)
m = np.concatenate(new_M)
f = LinearNDInterpolator(list(zip(x,y)), z)
f_mass = LinearNDInterpolator(list(zip(x,y)),m)
#%%

x = np.concatenate(new_T)
y = np.concatenate(new_R)
z = np.concatenate(new_A)
m = np.concatenate(new_M)
f = LinearNDInterpolator(list(zip(x,y)), z)
f_mass = LinearNDInterpolator(list(zip(x,y)),m)
#%%
num_iter = 5000
age_MC = []
mass_MC = []
newvals = []

for b in range(num_iter):
    
    R_err = np.random.normal(R_star,dR)
    T_err = np.random.normal(t_star,dt_star)
    MC_a = f(T_err, R_err)
    MC_mass = f_mass(T_err, R_err)
    
    # print("Age:",(10**MC_a)/1e6)
    # print("Mass:",MC_mass)
    newvals.append((T_err,R_err))    
    
    age_MC.append(MC_a)
    mass_MC.append(MC_mass)
    
    print("Iteration ", b+1, "out of ", num_iter)
#%%
#Monte Carlo simulation
num_iter = 5000
age_MC = []
mass_MC = []
newvals = []

for b in range(num_iter):
    
    L_err = np.random.normal(L_star,dL)
    T_err = np.random.normal(Teff_star,dT)
    MC_a = f(T_err, L_err)
    MC_mass = f_mass(T_err, L_err)
    
    # print("Age:",(10**MC_a)/1e6)
    # print("Mass:",MC_mass)
    newvals.append((T_err,L_err))    
    
    age_MC.append(MC_a)
    mass_MC.append(MC_mass)
    
    print("Iteration ", b+1, "out of ", num_iter)
    
#%%
new_temps = []
age_yr = []
rad_nl = []
mass_nl = []
for i in range(num_iter):
    t = newvals[i][0]
    a = (10**age_MC[i])/1000000
    r = newvals[i][1]
    mass = mass_MC[i]/(1.989e33)
    new_temps.append(t)
    age_yr.append(a)
    rad_nl.append(r)
    mass_nl.append(mass)
#%%
#converts data out of log space
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
MC_data = pd.read_csv("PARSEC_MC_906.csv", header = None)
lum_nl = MC_data[0]
new_temps = MC_data[1]
age_yr = MC_data[2]
mass_MC = MC_data[3]

#%%

samples = np.vstack([age_yr, new_temps, mass_nl, rad_nl])
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'

figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),range = [(18,36), (7300,7600),(1.3,1.8),(1.4,1.5)], labels = [r"$\rm Age [Myr]$", r"$T_{\rm eff} \rm [K]$",r"$Mass [M_{\odot}]$", r"$R_{\star} [\rm R_{\odot}]$"], show_titles = False,title_kwargs={"fontsize": 15}, max_n_ticks = 5 )
axes = np.array(figure.axes).reshape((4,4))

ageq = corner.quantile(samples[0], q = (0.16, 0.5,0.84))
massq = corner.quantile(samples[2], q = (0.16, 0.5,0.84))
teffq = corner.quantile(samples[1], q = (0.16, 0.5,0.84))
radq = corner.quantile(samples[3], q = (0.16, 0.5,0.84))

ageup = ageq[2]-ageq[1]
agedown = ageq[1]-ageq[0]
teffup = teffq[2]-teffq[1]
teffdown = teffq[1]-teffq[0]
radup = radq[2]-radq[1]
raddown = radq[1]-radq[0]
massup = massq[2]-massq[1]
massdown = massq[1]-massq[0]
axes[0, 0].set_title(r"$\rm Age~[Myr] =  %.1f ^{+%.1f}_{- %.1f}$" % (ageq[1], ageup, agedown), fontsize = 16, pad=10) 
axes[3, 3].set_title(r"$\rm Radius~[R_{\odot}] =  %.3f ^{+%.3f}_{- %.3f}$" % (radq[1], radup, raddown),fontsize = 16, pad = 10)
axes[1, 1].set_title(r"$T_{\rm eff}~[K] =  %d ^{+%d}_{- %d}$" % (teffq[1], teffup, teffdown),fontsize = 16, pad = 10)
axes[2, 2].set_title(r"$M_{\star}~[\rm M_{\odot}] =  %.2f ^{+%.2f}_{- %.2f}$" % (massq[1], massup, massdown),fontsize = 16, pad = 10)
#%%
#plots the corner plots
samples = np.vstack([age_yr, mass_MC, new_temps, lum_nl])
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['axes.autolimit_mode'] = 'round_numbers'

figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84),range = [(18,30),(1.530,1.565), (7300,7600),(5.4,6.0)],labels = [r"$\rm Age [Myr]$", r"$\rm Mass [ M_{\odot}$]", r"$T_{\rm eff} \rm [K]$", r"$L_{\star} [\rm L_{\odot}]$"], show_titles = False,title_kwargs={"fontsize": 15}, max_n_ticks = 5 )
axes = np.array(figure.axes).reshape((4,4))
# figure.xaxis.set_major_locator(ticker.MultipleLocator(base=2))  # x-axis ticks at multiples of 2
# figure.yaxis.set_major_locator(ticker.MultipleLocator(base=0.5))  # y-axis ticks at multiples of 0.5

ageq = corner.quantile(samples[0], q = (0.16, 0.5,0.84))
massq = corner.quantile(samples[1], q = (0.16, 0.5,0.84))
teffq = corner.quantile(samples[2], q = (0.16, 0.5,0.84))
lumq = corner.quantile(samples[3], q = (0.16, 0.5,0.84))

ageup = ageq[2]-ageq[1]
agedown = ageq[1]-ageq[0]
massup = massq[2]-massq[1]
massdown = massq[1]-massq[0]
teffup = teffq[2]-teffq[1]
teffdown = teffq[1]-teffq[0]
lumup = lumq[2]-lumq[1]
lumdown = lumq[1]-lumq[0]


axes[3, 0].set_xticks([20.0,22.5,25.0,27.5])
axes[3, 1].set_xticks([1.53,1.53,1.54,1.55,1.56])
axes[3, 2].set_xticks([7300,7350,7400,7450,7500,7550])
axes[3, 3].set_xticks([5.5,5.6,5.7,5.8,5.9])
axes[2, 0].xaxis.set_tick_params(labelbottom=False)
axes[1, 0].set_xticks([20.0,22.5,25.0,27.5])
axes[1, 0].xaxis.set_tick_params(labelbottom=False)
axes[0, 0].set_xticks([20.0,22.5,25.0,27.5])
axes[0, 0].xaxis.set_tick_params(labelbottom=False)
axes[2, 1].set_xticks([1.53,1.53,1.54,1.55,1.56])
axes[2, 1].xaxis.set_tick_params(labelbottom=False)
axes[1, 1].set_xticks([1.53,1.53,1.54,1.55,1.56])
axes[1, 1].xaxis.set_tick_params(labelbottom=False)
axes[2, 2].set_xticks([7300,7350,7400,7450,7500,7550])
axes[2, 2].xaxis.set_tick_params(labelbottom=False)

axes[1, 0].set_yticks([1.53,1.53,1.54,1.55,1.56])
axes[2, 0].set_yticks([7300,7350,7400,7450,7500,7550])
axes[3, 0].set_yticks([5.5,5.6,5.7,5.8,5.9])
axes[2, 1].set_yticks([7300,7350,7400,7450,7500,7550])
axes[2, 1].yaxis.set_tick_params(labelleft=False)
axes[3, 1].set_yticks([5.5,5.6,5.7,5.8,5.9])
axes[3, 1].yaxis.set_tick_params(labelleft=False)
axes[3, 2].set_yticks([5.5,5.6,5.7,5.8,5.9])
axes[3, 2].yaxis.set_tick_params(labelleft=False)




axes[0, 0].set_title(r"$\rm Age~[Myr] =  %.1f ^{+%.1f}_{- %.1f}$" % (ageq[1], ageup, agedown), fontsize = 16, pad=10) 
axes[1, 1].set_title(r"$\rm Mass~[M_{\odot}] =  %.3f ^{+%.3f}_{- %.3f}$" % (massq[1], massup, massdown),fontsize = 16, pad = 10)
axes[2, 2].set_title(r"$T_{\rm eff}~[K] =  %d ^{+%d}_{- %d}$" % (teffq[1], teffup, teffdown),fontsize = 16, pad = 10)
axes[3, 3].set_title(r"$L_{\star}~[\rm L_{\odot}] =  %.2f ^{+%.2f}_{- %.2f}$" % (lumq[1], lumup, lumdown),fontsize = 16, pad = 10)
# axes[3,0].xaxis.set_major_locator(ticker.MultipleLocator(base = 0.05))
# figure.savefig('PARSEC_corner_0627_V2.jpg')
figure.savefig('PARSEC_corner_1011.pdf')
#%%
#write the monte carlo values to a file to save
from numpy import savetxt
MC_L = np.array(lum_nl)
MC_T = np.array(new_temps)
MC_age = np.array(age_yr)
MC_mass = np.array(mass_MC)

dat_MC = np.array([MC_L, MC_T, MC_age, MC_mass])
dat_MC = dat_MC.T
savetxt('PARSEC_MC_906.csv', dat_MC, delimiter = ',')