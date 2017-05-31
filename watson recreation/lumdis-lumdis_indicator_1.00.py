from scipy import * 
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import LambdaCDM
from astropy.io import ascii
import astropy.units as u




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

def breakuptable(file, rows= [1,10000], columns=[0,10000]):
    raw = ascii.read(file)
    table = {}
    
    names = list(raw.colnames)
    #changes rows if specified
    rstart = rows[0]
    rend = rows[1]
    
    #changes the columns that will be added to the list if specified
    cstart = columns[0]
    cend = columns[1]
    names = names[cstart:cend] 
    
    for i in range(len(names)):
        current_column = list(raw[[names[i]][0]])
        current_name = current_column[0]
        table[current_name] = current_column[rstart:rend]
        
    return table
    

dic = breakuptable('CSVResults z.csv')
 
names = dic.keys()






######################################### ANALYZING DATA ####################################################



R = []
F = []
Z = []
lumdis_indicator = []

def lumdis(z, key):
    Ho = cosmology[key].H
    mo = cosmology[key].Om
    l = cosmology[key].Ok
    cos = LambdaCDM(H0=Ho, Om0=mo, Ode0=l)    
    return cos.luminosity_distance(z)




def getZ():
    for p in dic['z']:
        Z.append(float(p))
       


def getF(): #getting flux from luminosity 
    n = 0 
    for i in dic['log_L_agn']:
        z = Z[n]
        l = 10**float(i)
        dis = float((lumdis(z,'conversion'))/u.Mpc)
        f = (l)/(4 * pi * (dis ** 2))
        F.append(f)
        n += 1


def getR():
    c = 29979245800
    for x in dic['t_cent']:
        y = float(x)
        z = (y * c)
        R.append(z)
# broad line region radius: R = c(tau)


def getlumdis_indicator():
    n = 0
    lambd = 4861 #angstoms 
    for t in dic['t_cent']:
        t = float(t)
        indc = t/((F[n]*lambd)**.5)
        lumdis_indicator.append(indc)





getZ()
getR()
getF()

getlumdis_indicator()



######################################### MAKING GRAPH ######################################################




z = np.arange(0,4,.001)
norm = lumdis(z, "s")
bs0 = lumdis(z, "bs0")
bs1 = lumdis(z, "bs1")
bs2 = lumdis(z, "bs2")

fig, ax1 = plt.subplots()


t = np.array(lumdis_indicator)


        
ax1.loglog(z, norm, 'k', label="Omega m = .3089 Omega lambda = .6911 Omega k = 0 (Current Best Cosmolgy)")
ax1.loglog(z, bs0, 'r:', label="Omega m = .1 Omega lambda = .9 Omega k = 0")
ax1.loglog(z, bs1, 'b--', label="Omega m = 0 Omega lambda = 0 Omega k = 1")
ax1.loglog(z, bs2, 'g-.', label="Omega m = .5 Omega lambda = 1 Omega k = -.5")
ax1.set_xlabel("Z (Redshift)", fontsize = 14.5)
plt.title("Luminosity Distance vs Z For Different Cosmologies", fontsize = 20)
ax1.set_ylabel("Luminosity Distance (MParsecs)", fontsize = 14.5)

ax2 = ax1.twinx()
ax2.loglog(Z, t, 'r.')
ax2.set_ylabel("tau/sqrt(F)", fontsize = 14.5)

   

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




legend = plt.legend(loc="upper left")

"""
for label in legend.get_texts():
    label.set_fontsize('small')
"""

plt.show()


# Calculated luminosity distances are far too high
# Flux is wrong, can't use R should be distance to observer, not R of broad line region
