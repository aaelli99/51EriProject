import os
os.chdir("C:\\Users\\oxfor\\OneDrive\\Documents\\Ashley's Stuff\\School stuff\\Grad School\\Research")
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import pandas as pd

#%%
# f = open('match_Z0.017Y0.279OUTA1.74_F7_M001.600.txt','r')

data = pd.read_csv("match_Z0.017Y0.279OUTA1.74_F7_M001.500.txt", sep='\s+')

age = data['#']
teff = data['Mass']
mbol = data['logTe']

plt.plot(teff,mbol)
plt.axis([4.0,3.5,5,-5])
plt.xlabel('Teff')
plt.ylabel('Mag')


maxt = max(teff)
idx_maxt = np.argmax(teff)
print("Max age:", 10**(age[idx_maxt]))
