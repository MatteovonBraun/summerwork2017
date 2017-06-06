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

trim['pct_err_t'] = trim.apply(lambda row: percenterror(row), axis = 1)

idx =  trim.groupby(['name1'])['pct_err_t'].transform(min) == trim['pct_err_t']
full = trim[idx]

#################################################ANALISIS#######################################################

#using the 'clean2+ExtCorr' values on page 22 of Bentz et al. 2013
alpha = .549
beta =  1.559   #'K' in paper

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
    dis = float((lumdis(z, 'conversion'))/u.Mpc)*(3.08567758*10**24) #conversion to cm
    F = row['L_agn']/(4*m.pi*(dis**2))
    return F

full['F'] = full.apply(lambda row: flux(row), axis = 1)

def radius(row):
    #converts to seconds, then multiplied by c in cm
    r = row['t_cent']*86400*29979245800 
    return float(r)

full['r_cm'] = full.apply(lambda row: radius(row), axis = 1)

def radius_e(row):
    r_e = row['t_c_perr'] *86400*29979245800 
    return float(r_e)

full['r_cm_e'] = full.apply(lambda row: radius_e(row), axis = 1)

print full

def dis(row):
    c = 1/(m.sqrt(4*m.pi*row['F'])) #'c' stands for constant, not speed of light. 
    x = m.sqrt((row['r_cm']/( 10**beta ))**(1/alpha))  
    dis = c*x
    return dis

full['dis_mpc'] = full.apply(lambda row: dis(row), axis = 1)

#error in just Flux and tau
def dis_Ft_err(row): 
    err_F = 10**(row['log_L_agn_err'])
    err_r = row['r_cm_e']
    F = row['F']
    r = row['r_cm']
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+F
    ddF =  -1*m.sqrt((10**(-1*beta) * r)**(1/alpha))/(4* m.sqrt(m.pi)*F**(3/2))
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+t
    ddr = m.sqrt((10**(-1*beta)* r)**(1/alpha))/(4* m.sqrt(m.pi) * alpha * m.sqrt(F) * r)
    
    dis_err = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    return dis_err
              
full['dis_Ft_err'] = full.apply(lambda row: dis_Ft_err(row), axis = 1)

#error in flux, tau, alpha, and beta
def dis_Ftab_err(row):
    err_a = .0275
    err_b = .024
    F = row['F']
    r = row['r_cm']
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+a
    dda = -1*(m.sqrt((10**(-1*beta) * r)**(1/alpha)) * m.log(10**(-1*beta) * r))/(4 * m.sqrt(m.pi) * alpha**2 * m.sqrt(F))
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+b
    ddb =  -1*(m.log(10) * m.sqrt((10**(-1*beta) * r)**(1/alpha)))/(4 * m.sqrt(m.pi) * alpha * m.sqrt(F))
    
    dis_err_full = m.sqrt( row['dis_Ft_err']**2 + (dda*err_a)**2 + (ddb*err_b)**2 )
    return dis_err_full 

full['dis_Ftab_err'] = full.apply(lambda row: dis_Ftab_err(row), axis = 1)
    

######################################### MAKING GRAPH ######################################################

#Labels:

title = 'Luminosity  Distance  vs  Z  For  Different  Cosmologies'
axis_x =  "$Redshift$"
axis_y1 = 'Luminosity Distance $(MPc)$'
axis_y2 = r"$ \frac{\tau}{\sqrt{\lambda  F}}$"



z = np.arange(0,4,.001)
norm = lumdis(z, "s")
bs0 = lumdis(z, "bs0")
bs1 = lumdis(z, "bs1")
bs2 = lumdis(z, "bs2")

fig, ax1 = plt.subplots()

Z = full['z']
D = full['dis_mpc']
D_e= full['dis_Ft_err']
D_ee=full['dis_Ftab_err']
#t_ne=np.array(lumdis_indicator_ne) 
        
ax1.loglog(z, norm, 'k', label=r"$\Omega_m = .3089$  $\Omega_\lambda = .6911$  $\Omega_k = 0$  (Current Best Cosmolgy)", linewidth = 1)
ax1.loglog(z, bs0, 'r:', label=r"$\Omega_m = .1$  $\Omega_\lambda = .9 $  $\Omega_k = 0$", linewidth = 1)
ax1.loglog(z, bs1, 'y--', label=r"$\Omega_m = 0$   $\Omega_\lambda = 0$    $\Omega_k = 1$", linewidth = 1)
ax1.loglog(z, bs2, 'g-.', label=r"$\Omega_m = .5$  $\Omega_\lambda = 1$  $\Omega_k = -.5$", linewidth = 1)
ax1.set_xlabel(axis_x, fontsize = 14.5)
plt.title(title, fontsize = 15)
ax1.set_ylabel(axis_y1, fontsize = 14.5)

ax1.errorbar(Z, D,xerr=0, yerr = D_e, fmt = 'k.', ecolor='r'  , elinewidth= 1.4, label=r'Redshift Independant Distance $without$ $\alpha$ and $\beta$ error',zorder=10)
ax1.errorbar(Z, D,xerr=0, yerr = D_ee, fmt = 'k.', ecolor='#FF9900'  , elinewidth= 1, label=r'Redshift Independant Distance $with$ $\alpha$ and $\beta$ error',zorder=9)
#ax2.set_xscale('log', nonposx='clip')
#ax2.set_yscale('log', nonposy='clip')
#ax2.set_ylabel(axis_y2, fontsize = 19.5)

#ax1.axis([10**-3,4,4,1.1*10**5])

"""
In case you want to ha

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
