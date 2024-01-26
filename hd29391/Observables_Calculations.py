#Ashley Elliott
#Measured Observables for HD29391
import numpy as np
import matplotlib.pyplot as plt
import math as m
from zero_point import zpt
import pandas as pd

#%%
#Functions to easily calculate everything

#function to calculate linear radius from angular diameter
def ang_to_lin(theta,dtheta,D,dD):
    Rs = 6.957e8;               #meters in a solar radii
    pc_conv = 3.091e16          #meters in a parsec
    mas_conv = 206265*1000       #radians to mas conversion
    
    R = (theta*(D*pc_conv))/(2*mas_conv*Rs)
    # dR = R*np.sqrt((dtheta/theta)**2 + (dD/D)**2)                         #my equation same as Jonas just different derivation
    dR = (pc_conv/ (2*mas_conv*Rs)) * np.sqrt((dtheta*D)**2+(dD*theta)**2)  #jonas equation matches the first one
    return(print("Linear Radius: ", R ,"+/-", dR, "[R_solar]"))
    # return((R,dR))
#function to calculate temperature from bolemetric flux and angular diameter
def temp(Fbol,dF,theta,dtheta):
    T = 2341*(Fbol/theta**2)**(1/4)
    dT = np.sqrt(
          (2341 / (4 * Fbol**(3/4) * theta**(1/2)) * dF)**2 
        + (2341 * Fbol**(1/4) / (2* theta**(3/2)) * dtheta)**2
    )
    # dT = T*np.sqrt(((dF/Fbol)*(1/4))**2 + ((dtheta/theta)*(1/2))**2)
    return(print("Temperature: ", T ,"+/-",dT,"[K]"))
    # return((T,dT))
def luminosity(D,dD, Fbol, dF):
    pc_conv = 3.091e16          #meters in a parsec
    d = D*pc_conv*100;
    L = (4*np.pi*(d**2)*Fbol)/(3.846e33)
    dL = L*np.sqrt(((2*dD)/D)**2 + (dF/Fbol)**2)
    return(print("Luminosity: ", L ,"+/-",dL, "[L_solar]"))
    # return((L,dL))

def monte_D(p, p_err, num_it):
    parallaxes = np.random.normal(p, p_err, num_it)
    
    dist = 1/parallaxes
    avg = np.mean(dist)
    new_err = np.std(dist)
    
    # return(print(f"{avg:.4f} p/m {new_err:.4f}"))
    return((avg, new_err))
    
def monte(val, err, num_iter):
    rand_val = np.random.normal(val, err, num_iter)
    
    avg = np.mean(rand_val)
    new_err = np.std(rand_val)
    # return(print(f"{avg:.4f} p/m {new_err:.4f}"))
    return((avg, new_err))

def dist_calc(p, p_err, zpc):
    pc = (p+zpc)/1000
    D = 1/(pc)
    dD = D*((p_err/1000)/pc)
    print("Corrected parallax:", pc*1000)
    print("Distance:",D,"+/-", dD, "[pc]")
    return((D,dD))
#%%
zpt.load_tables()

# data = pd.read_csv('HD29391GAIAvals.csv')
phot_g_mean_mag = 5.141171
nu_eff_used_in_astrometry = 1.659
pseudocolour = 7.2
ecl_lat = -24.3053334255
astrometric_params_solved = 31
correction = zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)
#%%
p = 33.439                  #units in mas
pc_err = 0.0777             #units in mas
theta = 0.4497730           #angular diameter in mas
Fbol = 20.46            #bolometric flux in 10e-8 ergs/cm^2/s
dtheta = 0.004          #error in ang diam in mas
dF = 0.33               #bolometric flux error in ergs/cm^2/s
# MC_D = monte_D(pc/1000,pc_err/1000,10000)
# D = MC_D[0]
# dD = MC_D[1]
distance = dist_calc(p,pc_err,correction)
D = distance[0]
dD = distance[1]
lD = ang_to_lin(theta,dtheta,D,dD)
# MC_R = monte(lD[0], lD[1], 10000)
T = temp(Fbol, dF, theta, dtheta)
# MC_T = monte(T[0], T[1], 10000)
L = luminosity(D, dD, Fbol*(1e-8), dF*(1e-8))
# MC_L = monte(L[0], L[1], 10000)



# dD = 0.07               #error in distance in parsecs
# lD = ang_to_lin(theta,dtheta,D,dD)
# T = temp(Fbol,dF,theta, dtheta)
# L = luminosity(D,dD, Fbol*(1e-8), dF*(1e-8))


