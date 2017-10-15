from scipy import * 
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import LambdaCDM
from astropy.io import ascii
import astropy.units as u
import math as m 



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


print names



######################################### ANALYZING DATA ####################################################



R = []
F = []
F_e = []
Z = []
lumdis_indicator = []
lumdis_indicator_pe = []
lumdis_indicator_ne = []


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
    for l in dic['log_L_agn']:
        z = Z[n]
        l = float(l)
        
        L = 10**l
        L_e = (10**l)*(m.log(10))*(float(dic['log_L_agn_err'][n]))
        
        dis = float((lumdis(z,'conversion'))/u.Mpc)
        f = (L)/(4 * m.pi * (dis ** 2))
        f_e = (1/(4*m.pi*dis**2))*(L_e)
        
        
        F.append(f)
        F_e.append(f_e)
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
    lambd = 4861*(10**-7) #angstoms 
    for t in dic['t_cent']:
        t = float(t)
        t_pe = float(dic['t_c_perr'][n])
        #t_ne = float(dic['t_c_nerr'][n])
        
        f = float(F[n])
        f_e = float(F_e[n])
        
        
        
        indc = t/((f*lambd)**.5)
        indc_pe = m.sqrt(abs(((( 1/(m.sqrt(lambd*f)) )**2)*((t_pe)**2)) + ((((-1*lambd*t)/(2*(lambd*f)**(3/2))**2)*((f_e)**2)))))
        #indc_pe = m.sqrt(abs(((( 1/(m.sqrt(lambd*f)) )**2)*((t_pn)**2)) + ((((-1*lambd*t)/(2*(lambd*f)**(3/2))**2)*((f_e)**2)))))
         
        lumdis_indicator.append(indc)
        lumdis_indicator_pe.append(indc_pe)
        #lumdis_indicator_ne.append(indc_ne)



getZ()
getF()
getR()
getlumdis_indicator()


######################################### MAKING GRAPH ######################################################

#Labels:

title = 'Luminosity  Distance  vs  Z  For  Different  Cosmologies'
axis_x =  'Luminosity Distance $(MParsecs)$'
axis_y1 = "Z $(Redshift)$"
axis_y2 = r"$ \frac{\tau}{\sqrt{\lambda  F}}$"





z = np.arange(0,4,.001)
norm = lumdis(z, "s")
bs0 = lumdis(z, "bs0")
bs1 = lumdis(z, "bs1")
bs2 = lumdis(z, "bs2")

fig, ax1 = plt.subplots()


t = np.array(lumdis_indicator)
t_pe_f=np.array(lumdis_indicator_pe)
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
