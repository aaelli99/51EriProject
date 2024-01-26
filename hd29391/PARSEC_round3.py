#Another rewrite of code
#Ashley Elliott
import os

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
#66_88 is the old set of data
#f = open("PARSEC_data_66_88_new.txt", 'r')
#f_new = open("PARSEC_data_edit_66_88_new.txt", 'a')

iso_dir = os.path.join('hd29391', 'Isochrone Data')
old_data = os.path.join('hd29391', 'old_data')
parsec_data = os.path.join(old_data, 'parsecdatafiles')

f1 = open(os.path.join(parsec_data, 'PARSEC_data_7_4.txt'),'r')
f_new1 = open(os.path.join(parsec_data, "PARSEC_data_edit_7_4.txt"),'a')

f_new1.write("Zini MH logAge Mini int_IMF Mass logL logTe logg label McoreTP C_O period0  period1  period2  period3  period4  pmode  Mloss  tau1m   X   Y   Xc  Xn  Xo  Cexcess  Z 	 mbolmag  Umag    Bmag    Vmag    Rmag    Imag    Jmag    Hmag    Kmag \n")
no = "#"
for line in f1:
    if not no in line:
        f_new1.write(line)
        f_new1.write("\n")

f_new1.close()      
f1.close()
#%%
#Cell reads in the data and determines the beginning and ending indices for each age in order to sort the data into separate arrays based on the age it goes with

#data=pd.read_csv("PARSEC_data_edit_66_88_new.txt", sep='\s+')
data = pd.read_csv("PARSEC_data_edit_0507.txt", sep='\s+')

#removing columns I wont need
data.drop(['label','McoreTP','C_O','period0','period1', 'period2','period3','period4','pmode','Mloss', 'tau1m','X','Y','Xc','Xn','Xo','Cexcess','Z','mbolmag','Umag', 'Bmag','Vmag','Rmag','Imag','Jmag','Hmag'],axis=1, inplace=True)

logL = []
logTeff = []
mass = []
#age = [6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6]   #ages for the 66_88 data set
age = np.round(np.arange(6.5,7.95,0.05),2)
#age = 6.6          #single isochrone age
lens = np.empty(len(age))

for ii in range(len(age)):
    idx = 0
    for i in range(len(data['logL'])):
        if np.round(data["logAge"][i],2) == age[ii]:
            # print("Data is: ", data["logAge"][i], "with index ", i)
            # print("Age is: ", age[ii])
            idx = idx+1
            #print(idx)
        lens[ii] = idx
    
    beg_pos = np.empty(len(age))
    end_pos = np.empty(len(age))
    #print(lens[ii])
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
    logL.append(lums)
    logTeff.append(temps)
    mass.append(Mass)
    


#%%

#%%
data = np.genfromtxt("PARSEC_data_66_88_new.txt", names=True)
df = pd.DataFrame(data)

df_dict = {age: df[df["logAge"]==age] for age in df["logAge"].unique()}

final_df = pd.DataFrame()

vals = ["logL", "logTe", "logg", "MH"]
for age, frame in df_dict.items():
    index = pd.MultiIndex.from_product(([age], frame.Mass), names=['age', 'mass'])
    new_df = pd.DataFrame(index=index)
    for val in vals:
        new_df[val] = frame[val].to_numpy()
        
    final_df = final_df.append(new_df)

