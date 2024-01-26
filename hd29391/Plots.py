#Ashley Elliott
#Plotting V^2 for HD29391
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.special as ss
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from astropy.io import ascii,fits
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
from scipy.stats import binned_statistic
from lmfit import Model
plt.rcParams['text.usetex'] = True
import pandas as pd
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

def v2(spat_freq, theta):
    mu = 0.4516
    theta = theta/(206265*1000)
    x = np.pi*spat_freq*theta
    return ((((1-mu)/2)+(mu/3))**(-1)*((1-mu)*(ss.jv(1,x)/x)+mu*np.sqrt(np.pi/2)*(ss.jv(3/2,x)/x**(3/2))))**2
#%% reading in the data and calculating the model
# data1 =ascii.read('HD29391_LDCFit_Redo_0623.txt')
data1 =ascii.read('HD29391_LDDFit_CorrectLDC.txt')
spat_freq = data1['col1']
sf = np.linspace(0.000000001,6e8,10000)
vis2 = data1['col2']
d_vis2 = data1['col3']

mu = 0.4516
# theta = 0.45078965
theta = 0.4497730
model = V2(sf,theta, mu)

mod_vals = np.array([sf,model])
mod_vals = mod_vals.T
real_vals = np.array([spat_freq,vis2])
real_vals = real_vals.T


#%%
new = sorted((i,j,k) for i,j,k in zip(spat_freq, vis2, d_vis2))
real_vals_sort = np.array(new)

df = pd.DataFrame({'Spat Freq':real_vals_sort[:,0], 'Vis' : real_vals_sort[:,1], 'Yerr': real_vals_sort[:,2]})
mu_p = 0.4516
mu_p_arr = np.ones(len(df['Vis']))*mu_p
#%%
new = sorted((i,j,k) for i,j,k in zip(spat_freq, vis2, d_vis2))
real_vals_sort = np.array(new)

min_sf_val = real_vals_sort[:,0].min()
max_sf_val = real_vals_sort[:,0].max()
num_bins = round((max_sf_val-min_sf_val)/(10**7))

bins = np.histogram_bin_edges(real_vals_sort[:,0], int(num_bins))
d = np.digitize(real_vals_sort[:,0],bins)
avg_sf = np.zeros([int(num_bins)+1])
avg_v2 = np.zeros([int(num_bins)+1])
avg_yerr = np.zeros([int(num_bins)+1])
count = np.zeros([int(num_bins)+1])

for i in range(int(num_bins)+1):
    tot_sf = 0
    tot_v2 = 0
    tot_w = 0
    for j in range(len(real_vals_sort[:,0])):
        if (d[j] == i+1):
            tot_sf = tot_sf + real_vals_sort[j,0]
            weight = 1/((real_vals_sort[j,2])**2)
            tot_v2 = tot_v2 + weight*real_vals_sort[j,1]
            tot_w = tot_w + weight
            count[i] = count[i] + 1
            
    avg_sf[i] = tot_sf/count[i]
    # print("Average Spat Freq:", avg_sf[i])
    avg_v2[i] = tot_v2/tot_w
    # print("Average Vis:", avg_v2[i])
    avg_yerr[i] = 1/(np.sqrt(tot_w))
    # print("Total Error:", tot_w)
    # print("Average Error:", avg_yerr[i])
    # print("Counts:", count[i])

#%% Calculating the residuals
avg_model = V2(avg_sf,theta,mu)
model_n = V2(spat_freq,theta,mu)

# avg_res = (abs(avg_v2-avg_model)/(avg_v2))*100
# res = (abs(vis2-model_n)/ (vis2))*100
avg_res = (avg_v2-avg_model)
res = (vis2-model_n)
avg_res_err = avg_yerr/avg_v2
res_err = real_vals_sort[:,2]/real_vals_sort[:,1]

#%%
#Gail and M. Simon's work

v2_gs = np.array([0.849,0.819,0.860,0.769,0.776,0.766,0.775,0.773,0.719,0.791,0.729])
v2_gs_err = np.array([0.013,0.012,0.014,0.021,0.022,0.021,0.021,0.022,0.020,0.022,0.020])
B_gs = np.array([313.26,313.44,312.45,311.41,313.41,312.20,310.15,311.11,313.43,313.31,311.69])    #baseline in m
lam = np.array([2.133,2.133,2.133,1.673,1.673,1.673,1.673,1.673,1.673,1.673,1.673])*(1e-6)   #wavelength in m
sf_gs = (B_gs/lam)

#%%
v2_c = np.array([0.710,0.7553,0.640])
v2_c_err = np.array([0.036,0.0338,0.027])
sf_c = np.array([1.838,1.888,1.917])*(1e8)
lam_c = np.array([1.673,1.673,1.673])*(1e-6)
B_c = lam_c*sf_c
# sf_c = (B_c/lam_c)

# model_c = V2(sf_c, theta, mu)
# res_c = (v2_c-model_c)
# res_c_err = v2_c_err/v2_c
mu_c = 0.2582
mu_c_arr = np.ones(len(v2_c))*mu_c
popt, pcov = curve_fit(v2, sf_c, v2_c, bounds = (0.00000001, np.inf), sigma = v2_c_err, absolute_sigma = False)
SE = np.sqrt(np.diag(pcov))
spatfreq = np.linspace(0.000000001,6e8,10000)
fit_v = v2(spatfreq, popt[0])
plt.plot(sf_c, v2_c,'.')
plt.plot(spatfreq,fit_v,'k')
#%%

full_v2 = np.concatenate((v2_c,df['Vis'].values))
full_v2_err = np.concatenate((v2_c_err,df['Yerr'].values))
full_sf = np.concatenate((sf_c, df['Spat Freq'].values))
full_mu = np.concatenate((mu_c_arr, mu_p_arr))
#%%
gmodel = Model(v2)
print(f'parameter names: {gmodel.param_names}')
print(f'independent variables: {gmodel.independent_vars}')

result = gmodel.fit(full_v2,spat_freq = full_sf, theta = 0.2)

print(result.fit_report())
#%%
params, covar = curve_fit(v2, full_sf, full_v2, bounds = (0.000001, np.inf),sigma = full_v2_err, absolute_sigma = False)
SE = np.sqrt(np.diag(covar))
fit_vis = v2(full_sf, params[0])

# # chi_sqr = np.sum(((full_v2-fit_vis)/(SE))**2)
# # DoF = len(full_v2)-1
# # chi_sqr_r = chi_sqr/DoF
plt.plot(full_sf, full_v2, '.')
# # plt.errorbar(full_sf, full_v2, yerr = full_v2_err, '.')
plt.plot(full_sf, fit_vis, 'k')
plt.show()
#%% plotting
plt.rcParams.update({'font.size': 12})
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
f, (a0,a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios': [10,3]}, sharex = True)
a0.plot(real_vals_sort[:,0],real_vals_sort[:,1],'.',color = 'gray', alpha = 0.2)
a0.plot(sf, model,'--')
a0.errorbar(real_vals_sort[:,0],real_vals_sort[:,1],yerr = real_vals_sort[:,2],linestyle = "None", linewidth = 0.5, color = 'gray', capsize = 3, alpha = 0.2)
a0.plot(avg_sf,avg_v2,'.k', label = 'PAVO')
a0.errorbar(avg_sf,avg_v2,yerr = avg_yerr, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
a0.plot(sf_c, v2_c, 's', markersize = 3, color = 'green', label = 'CLASSIC')
a0.errorbar(sf_c, v2_c, yerr = v2_c_err, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'green')
# a0.plot(sf_gs,(v2_gs), '.', color = 'blue', label = 'M. Simon and G.H. Schafer')
# a0.errorbar(sf_gs, v2_gs, yerr = v2_gs_err, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'blue')
a1.plot(spat_freq,res, '.', color = 'gray', alpha = 0.2)
a1.plot(avg_sf, avg_res,'.k')
a1.errorbar(spat_freq,res,yerr =real_vals_sort[:,2] ,linestyle = "None", linewidth = 0.5, color = 'gray', capsize = 3, alpha = 0.2)
a1.errorbar(avg_sf,avg_res,yerr = avg_yerr, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'k')
a1.plot(sf_c, res_c, 's', color = 'green', markersize = 3)
a1.errorbar(sf_c, res_c, yerr = res_c_err, linestyle = 'None', linewidth = 0.5, capsize = 3, color = 'green')
plt.subplots_adjust(wspace = 0, hspace = 0)
plt.xlabel(r'$\rm Spatial$ $\rm frequency$ [$\rm rad^{-1}$]')
a0.set_ylabel(r'$V^2$', labelpad = 10)
a1.set_ylabel(r'$\rm Residual$',labelpad = 3)
a1.axhline(y = 0)
# plt.yticks([-20,0,50,100])
plt.yticks([-.25,0,0.25])
plt.xlim([0,6e8])
a0.legend()
a0.tick_params(axis = 'x', labelbottom = False)
# eq = r"$\mu = 0.4516$"
eq1 = r"$\theta_{LD} = 0.450 \pm 0.004 \rm [mas]$"
# a0.text(2e7, 0.18, eq, {'color': 'k', 'fontsize': 13})
a0.text(2e7, 0.10, eq1, {'color': 'k', 'fontsize':13})
a0.xaxis.set_minor_locator(AutoMinorLocator())
a0.yaxis.set_minor_locator(AutoMinorLocator())
a1.xaxis.set_minor_locator(AutoMinorLocator())
a1.yaxis.set_minor_locator(AutoMinorLocator())
f.savefig('HD29391_plot_1215.pdf')

#%%
plt.plot(real_vals_sort[321,0], real_vals_sort[321,1],'.g')
plt.plot(avg_sf[29], avg_v2[29],'.b')
#%% Uniform Limb Darkening not needed for now
# data2 =ascii.read('HD29391_UDFit.txt')
# spat_freq_U = data2['col1']
# sf_U = np.linspace(0,6e8,1000000)
# vis2_U = data2['col2']
# d_vis2_U = data2['col3']

# mu_U = 0
# theta = 0.43233056
# model_U = V2(sf_U,theta, mu_U)
# plt.figure(2)
# plt.plot(spat_freq_U,vis2_U,'k.',sf_U,model_U)
# plt.errorbar(spat_freq_U,vis2_U,yerr = d_vis2_U,linestyle = "None", linewidth = 0.5, color = 'black', capsize = 3)
# plt.xlabel(r'Spatial frequency [rad^{-1}]')
# plt.ylabel(r'$V^2$')
# plt.title(r'Uniform Limb Darkening Coeeficient $\mu = 0$, $\theta_{UDD} = 0.432\pm 0.0011$ [mas]')
