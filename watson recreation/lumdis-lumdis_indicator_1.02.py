from scipy import * 
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import LambdaCDM
import astropy.units as u
import math as m 
import pandas as pd 



#############################################SETUP###################################################

class cosmo(object):
    def __init__(self, H, Ok, Om, Ode):
        self.H = H
        self.Ok = Ok
        self.Om = Om
        self.Ode = Ode
        



cosmology = {"s":cosmo(70, 0, .3089, .6911), 
             "bs0":cosmo(70, 0,.1,.9), 
             "bs1":cosmo(70, 1, 0, 0), 
             "bs2":cosmo(70, -.5,.5 , 1),
             "conversion":cosmo(71,0, .3,.7)  #values gotten from http://adsabs.harvard.edu/abs/2015PASP..127...67B
            } 


#the following is the pandas setup, see commit link for details. 
#https://stackoverflow.com/questions/15705630/python-getting-the-row-which-has-the-max-value-in-groups-using-groupby
raw = pd.read_csv('CSVResults z.csv')
raw.sort_values(by='name1')

trim = pd.DataFrame()
trim = raw[['name1','log_L_agn','log_L_agn_err','t_cent', 't_c_perr', 'z']]

def percenterror(row):
    return (row['t_c_perr']/row['t_cent'])

trim['pct_err'] = trim.apply(lambda row: percenterror(row), axis = 1)

idx =  trim.groupby(['name1'])['pct_err'].transform(min) == trim['pct_err']
full = trim[idx]

##########################################ANALISIS####################################################


def lumdis(z, key):
    Ho = cosmology[key].H
    mo = cosmology[key].Om
    l = cosmology[key].Ok
    cos = LambdaCDM(H0=Ho, Om0=mo, Ode0=l)    
    return cos.luminosity_distance(z)


#source: https://stackoverflow.com/questions/26886653/pandas-create-new-column-based-on-values-from-other-columns

def lum10(row):
    L = 10**row['log_L_agn']
    return L

full['L_agn'] = full.apply(lambda row: lum10(row), axis = 1)

def flux(row):
    z = row['z']
    dis = float((lumdis(z, 'conversion'))/u.Mpc) #Im just realising that the units here might be wrong. 
    F = row['L_agn']/(4*m.pi*(dis**2))
    return F

full['F'] = full.apply(lambda row: flux(row), axis = 1)


def indc(row):
    lambd = 4861 
    t = row['t_cent']
    F = row['F']
    indc = t/(m.sqrt(lambd*F)) 
    return indc 

full['indc'] = full.apply(lambda row: indc(row), axis = 1)
    

######################################### MAKING GRAPH ######################################################

#Labels:

title = 'Luminosity  Distance  vs  Z  For  Different  Cosmologies'
axis_x =  "Z $(Redshift)$"
axis_y1 = 'Luminosity Distance $(MParsecs)$'
axis_y2 = r"$ \frac{\tau}{\sqrt{\lambda  F}}$"



z = np.arange(0,4,.001)
norm = lumdis(z, "s")
bs0 = lumdis(z, "bs0")
bs1 = lumdis(z, "bs1")
bs2 = lumdis(z, "bs2")

fig, ax1 = plt.subplots()

Z = full['z']
t = full['indc']
t_pe_f= 0
#t_ne=np.array(lumdis_indicator_ne) 
        
ax1.loglog(z, norm, 'k', label="Omega m = .3089 Omega lambda = .6911 Omega k = 0 (Current Best Cosmolgy)")
ax1.loglog(z, bs0, 'r:', label="Omega m = .1 Omega lambda = .9 Omega k = 0")
ax1.loglog(z, bs1, 'b--', label="Omega m = 0 Omega lambda = 0 Omega k = 1")
ax1.loglog(z, bs2, 'g-.', label="Omega m = .5 Omega lambda = 1 Omega k = -.5")
ax1.set_xlabel(axis_x, fontsize = 14.5)
plt.title(title, fontsize = 15)
ax1.set_ylabel(axis_y1, fontsize = 14.5)

ax2 = ax1.twinx()
ax2.set_xscale('log', nonposx='clip')
ax2.set_yscale('log', nonposy='clip')
ax2.errorbar(Z, t,xerr=0, yerr = t_pe_f, fmt = 'r.', elinewidth= .3)
ax2.set_ylabel(axis_y2, fontsize = 19.5)

ax2.axis([.001,4,4*10**-22,9*10**-18]) #this lines up the data 
   

"""
USE THIS IN THE FUTURE TO DO THE LINEAR TRANSFORMATION

def plotLR():
    x = L
    y = R
    
    #makes the bounds general so the data looks nice. 
    xgap = (.05*(max(x)-min(x)))
    ygap = (.07*(max(y)-min(y)))
    
    xlow = min(x) - xgap
    xhigh = max(x) + xgap
    ylow = min(y) - ygap
    yhigh = max(y) +ygap
    
    
    plt.plot(L,R, 'g.')
    plt.axis([xlow,xhigh,ylow,yhigh])
    plt.xlabel('log(L)')
    plt.ylabel('log(R)')
    plt.title('Luminosity vs Radius')
    plt.show()
"""




legend = ax1.legend(loc="upper left")


for label in legend.get_texts():
    label.set_fontsize('small')


plt.show()


# Calculated luminosity distances are far too high
# Flux is wrong, can't use R should be distance to observer, not R of broad line region
