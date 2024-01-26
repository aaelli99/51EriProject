#Ashley Elliott
#Plotting V^2 for HD29391
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii,fits
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import pandas as pd
import scipy
import scipy.special as ss
plt.rcParams['text.usetex'] = True

#%% Visibility Squared Functions
def V2(spat_freq,theta,mu):
    #############################################################################################
    #    function calculates the visibility squared for a star given a baseline, observational  # 
    #    wavelength, and angular diameter                                                       #
    #    Variables:                                                                             #
    #    mu = linear limb darkening coefficient                                                 #
    #    lambda = observational wavelength                                                      #
    #    B = baseline                                                                           #
    #    theta = angular diameter in radians                                                    #
    #############################################################################################
    theta = theta/(206265*1000) 
    x = (np.pi*spat_freq*theta)
    v = (((1-mu)/2)+(mu/3))**(-1)*((1-mu)*(ss.jv(1,x)/x)+mu*np.sqrt(np.pi/2)*(ss.jv(3/2,x)/x**(3/2)))
    return v**2

#%%
sf = np.linspace(0.00000000,6e8,10000)
mu = 0.4516
theta = 0.45078965
model = V2(sf,theta, mu)

#%%
df = pd.read_csv("HD29391_LDDFit_Correct.csv")
spat_freq = df['Spatial Frequency'].to_numpy()
v2 = df['V^2'].to_numpy()
yerr = df['Delta V^2'].to_numpy()

bin_size = 10**7
bin_ids = spat_freq//bin_size

binned_sf = []
binned_v2 = []
binned_yerr = []
for bin_id in np.unique(bin_ids):
    binned_sf.append(np.average(spat_freq[bin_ids==bin_id]))
    binned_v2.append(np.average(v2[bin_ids==bin_id], weights = 1/(yerr[bin_ids==bin_id])**2))
    N = sum(bin_ids==bin_id)
    binned_yerr.append(1/np.sqrt(sum(1/yerr[bin_ids==bin_id]**2)))

#%%
model_n = V2(spat_freq,theta, mu)
avg_model = V2(np.array(binned_sf),theta, mu)

#fractional difference residiual
avg_res_FD = (binned_v2-avg_model)/avg_model
res_FD = (v2-model_n)/model_n

#O-C residual
avg_res_OC = binned_v2-avg_model
res_OC = v2-model_n

#Percentage residual
avg_res_P = ((binned_v2-avg_model)/avg_model)*100
res_P = ((v2-model_n)/model_n)*100

#Sigma residual
avg_res_S = (binned_v2-avg_model)/binned_yerr
res_S = (v2-model_n)/yerr


#%% Fractional Difference Plot
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
f, (a0,a1) = plt.subplots(2,1, sharex = True, gridspec_kw = {'height_ratios': [10,3]})
a0.plot(spat_freq,v2,'.',color = 'gray', alpha = 0.2)
a0.plot(sf, model,'--')
a0.errorbar(spat_freq,v2,yerr = yerr,linestyle = "None", linewidth = 0.5, color = 'gray', capsize = 3, alpha = 0.2)
a0.plot(binned_sf,binned_v2,'.k')
a0.errorbar(binned_sf,binned_v2,yerr = binned_yerr, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
a1.plot(spat_freq,res_FD, '.', color = 'gray', alpha = 0.2)
a1.plot(binned_sf, avg_res_FD,'.k')
a1.errorbar(binned_sf,avg_res_FD,yerr = (np.array(binned_yerr)/np.array(binned_v2)), linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
plt.subplots_adjust(wspace = 0, hspace = 0)
plt.xlabel(r'$\rm Spatial$ $\rm frequency$ [$\rm rad^{-1}$]')
a0.set_ylabel(r'$V^2$', labelpad = 10)
a1.set_ylabel(r'$\rm Residual$',labelpad = 3)
a1.axhline(y = 0)
#plt.yticks([-0.5, 0,0.5])
plt.xlim([0,6e8])
a0.tick_params(axis = 'x', labelbottom = False)
eq = r"$\mu = 0.4516$"
eq1 = r"$\theta_{LD} = 0.451 \pm 0.0012 \rm [mas]$"
a0.text(0.5e8, 0.18, eq, {'color': 'k', 'fontsize': 12})
a0.text(0.5e8, 0.10, eq1, {'color': 'k', 'fontsize':12})
a0.xaxis.set_minor_locator(AutoMinorLocator())
a0.yaxis.set_minor_locator(AutoMinorLocator())
a1.xaxis.set_minor_locator(AutoMinorLocator())
a1.yaxis.set_minor_locator(AutoMinorLocator())
f.savefig('HD29391_plot_FD_1.pdf')


#%% Percent Difference Plot
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
f, (a0,a1) = plt.subplots(2,1, sharex = True, gridspec_kw = {'height_ratios': [10,3]})
a0.plot(spat_freq,v2,'.',color = 'gray', alpha = 0.1)
a0.plot(sf, model,'--')
a0.errorbar(spat_freq,v2,yerr = yerr,linestyle = "None", linewidth = 0.5, color = 'gray', capsize = 3, alpha = 0.2)
a0.plot(binned_sf,binned_v2,'.k')
a0.errorbar(binned_sf,binned_v2,yerr = binned_yerr, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
a1.plot(spat_freq,res_P, '.', color = 'gray', alpha = 0.2)
a1.plot(binned_sf, avg_res_P,'.k')
a1.errorbar(binned_sf,avg_res_P,yerr = (np.array(binned_yerr)/np.array(binned_v2))*100, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
plt.subplots_adjust(wspace = 0, hspace = 0)
plt.xlabel(r'$\rm Spatial$ $\rm frequency$ [$\rm rad^{-1}$]')
a0.set_ylabel(r'$V^2$', labelpad = 10)
a1.set_ylabel(r'$\rm Residual [\%]$',labelpad = 3)
a1.axhline(y = 0)
#plt.yticks([-50, 0,50])
plt.xlim([0,6e8])
a0.tick_params(axis = 'x', labelbottom = False)
eq = r"$\mu = 0.4516$"
eq1 = r"$\theta_{LDD} = 0.451 \pm 0.0012 \rm [mas]$"
a0.text(0.5e8, 0.18, eq, {'color': 'k', 'fontsize': 12})
a0.text(0.5e8, 0.10, eq1, {'color': 'k', 'fontsize':12})
a0.xaxis.set_minor_locator(AutoMinorLocator())
a0.yaxis.set_minor_locator(AutoMinorLocator())
a1.xaxis.set_minor_locator(AutoMinorLocator())
a1.yaxis.set_minor_locator(AutoMinorLocator())
#f.savefig('HD29391_plot.eps', format='eps')



#%% Sigma Offset Plot
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
f, (a0,a1) = plt.subplots(2,1, sharex = True, gridspec_kw = {'height_ratios': [10,3]})
a0.plot(spat_freq,v2,'.',color = 'gray', alpha = 0.2)
a0.plot(sf, model,'--')
a0.errorbar(spat_freq,v2,yerr = yerr,linestyle = "None", linewidth = 0.5, color = 'gray', capsize = 3, alpha = 0.2)
a0.plot(binned_sf,binned_v2,'.k')
a0.errorbar(binned_sf,binned_v2,yerr = binned_yerr, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
a1.plot(spat_freq,res_S, '.', color = 'gray', alpha = 0.2)
#a1.plot(binned_sf, avg_res_S,'.k')
plt.subplots_adjust(wspace = 0, hspace = 0)
plt.xlabel(r'$\rm Spatial$ $\rm frequency$ [$\rm rad^{-1}$]')
a0.set_ylabel(r'$V^2$', labelpad = 10)
a1.set_ylabel(r'$\rm Sigma Offset$',labelpad = 3)
a1.axhline(y = 0)
a1.axhline(y = 1,linestyle = '-.')
a1.axhline(y = 3,linestyle = '--')
a1.axhline(y = -1,linestyle = '-.')
a1.axhline(y = -3,linestyle = '--')
#plt.yticks([-1, 0,5])
plt.xlim([0,6e8])
a0.tick_params(axis = 'x', labelbottom = False)
eq = r"$\mu = 0.4516$"
eq1 = r"$\theta_{LDD} = 0.451 \pm 0.0012 \rm [mas]$"
a0.text(0.5e8, 0.18, eq, {'color': 'k', 'fontsize': 12})
a0.text(0.5e8, 0.10, eq1, {'color': 'k', 'fontsize':12})
a0.xaxis.set_minor_locator(AutoMinorLocator())
a0.yaxis.set_minor_locator(AutoMinorLocator())
a1.xaxis.set_minor_locator(AutoMinorLocator())
a1.yaxis.set_minor_locator(AutoMinorLocator())
#f.savefig('HD29391_plot.eps', format='eps')


#%%Uniform Limb Darkening not needed for now
data2 =ascii.read('HD29391_UDFit.txt')
spat_freq_U = data2['col1']
sf_U = np.linspace(0,6e8,1000000)
vis2_U = data2['col2']
d_vis2_U = data2['col3']

mu_U = 0
theta = 0.43233056
model_U = V2(sf_U,theta, mu_U)
plt.figure(2)
plt.plot(spat_freq_U,vis2_U,'k.',sf_U,model_U)
plt.errorbar(spat_freq_U,vis2_U,yerr = d_vis2_U,linestyle = "None", linewidth = 0.5, color = 'black', capsize = 3)
plt.xlabel(r'Spatial frequency [rad^{-1}]')
plt.ylabel(r'$V^2$')
plt.title(r'Uniform Limb Darkening Coeeficient $\mu = 0$, $\theta_{UDD} = 0.432\pm 0.0011$ [mas]')
