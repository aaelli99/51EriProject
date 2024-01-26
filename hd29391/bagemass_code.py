#Ashley Elliott
#bagemass isochrone data for HD 29391
import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research\\bagemass stuff")
#import read_mist_models
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import corner
# import arviz as az

#%%
#data = pd.read_csv("chain_0208.dat",  sep='\s+')
data = pd.read_csv("chain_0906.dat", sep = '\s+')
mass = data['Age'].to_numpy()
age = data['Step'].to_numpy()
temp = data['DeltaR'].to_numpy()
lum = data['Teff'].to_numpy()
#samples = np.vstack([age,mass,temp])
samples = np.vstack([age,mass,temp, lum])

mean_age = np.mean(age)
mean_temp = np.mean(temp)
mean_mass = np.mean(mass)
#%%
#converts data out of log space

age_yr = age/(1e-3)
lum_nl = 10**(lum)

samples = np.vstack([age_yr,mass,temp, lum_nl])
#%%
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 20})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
figure = corner.corner(samples.T,quantiles=(0.16, 0.5,0.84), labels = [r"$\rm Age [Myr]$", r"$\rm Mass [ M_{\odot}$]", r"$T_{eff} \rm [K]$", r"$L_{\star} [\rm L_{\odot}]$"], show_titles = True, title_fmt = ".3f",title_kwargs={"fontsize": 15})
axes = np.array(figure.axes).reshape((4,4))

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

axes[3, 0].set_xticks([30,35,40,45,50,55])                  #age
axes[3, 1].set_xticks([1.56,1.57,1.58,1.59,1.60,1.61])      #mass
axes[3, 2].set_xticks([7250,7300,7350,7400,7450])           #teff
axes[3, 3].set_xticks([5.4,5.6,5.8,6.0,6.2,6.4])            #lum
axes[2, 0].set_xticks([30,35,40,45,50,55])
axes[2, 0].xaxis.set_tick_params(labelbottom=False)
axes[1, 0].set_xticks([30,35,40,45,50,55])
axes[1, 0].xaxis.set_tick_params(labelbottom=False)
axes[0, 0].set_xticks([30,35,40,45,50,55])
axes[0, 0].xaxis.set_tick_params(labelbottom=False)
axes[2, 1].set_xticks([1.56,1.57,1.58,1.59,1.60,1.61])
axes[2, 1].xaxis.set_tick_params(labelbottom=False)
axes[1, 1].set_xticks([1.56,1.57,1.58,1.59,1.60,1.61])
axes[1, 1].xaxis.set_tick_params(labelbottom=False)
axes[2, 2].set_xticks([7250,7300,7350,7400,7450])
axes[2, 2].xaxis.set_tick_params(labelbottom=False)

axes[1, 0].set_yticks([1.56,1.57,1.58,1.59,1.60,1.61])      #mass
axes[2, 0].set_yticks([7250,7300,7350,7400,7450])           #teff
axes[3, 0].set_yticks([5.4,5.6,5.8,6.0,6.2,6.4])            #lum
axes[2, 1].set_yticks([7250,7300,7350,7400,7450])
axes[2, 1].yaxis.set_tick_params(labelleft=False)
axes[3, 1].set_yticks([5.4,5.6,5.8,6.0,6.2,6.4])
axes[3, 1].yaxis.set_tick_params(labelleft=False)
axes[3, 2].set_yticks([5.4,5.6,5.8,6.0,6.2,6.4])
axes[3, 2].yaxis.set_tick_params(labelleft=False)

axes[0, 0].set_title(r"$\rm Age~[Myr] =  %d ^{+%d}_{- %d}$" % (ageq[1], ageup, agedown), fontsize = 16, pad=10) 
axes[1, 1].set_title(r"$\rm Mass~[M_{\odot}] =  %.2f ^{+%.2f}_{- %.2f}$" % (massq[1], massup, massdown),fontsize = 16, pad = 10)
axes[2, 2].set_title(r"$T_{\rm eff}~[K] =  %d ^{+%d}_{- %d}$" % (teffq[1], teffup, teffdown),fontsize = 16, pad = 10)
axes[3, 3].set_title(r"$L_{\star}~[\rm L_{\odot}] =  %.2f ^{+%.2f}_{- %.2f}$" % (lumq[1], lumup, lumdown),fontsize = 16, pad = 10)
figure.savefig('GARSTEC_corner_1011.pdf')
# value1 = np.mean(samples.T,axis=0)

# for i in range(3):
#     ax = axes[i,i]
#     ax.axvline(value1[i], color = "r")
    
# for yi in range(3):
#     for xi in range(yi):
#         ax1 = axes[yi, xi]
#         ax1.axvline(value1[xi], color="r")
#         ax1.axhline(value1[yi], color="r")
#         ax1.plot(value1[xi], value1[yi], "sr")

