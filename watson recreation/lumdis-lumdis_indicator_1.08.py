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

"""
raw1 = pd.read_csv('mbh_measurements.csv')

join = pd.DataFrame({'name1':raw1['name'],
                     
                    
        
                        })
"""


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

def log_flux(row):
    z = row['z']
    log_dis = m.log10(float((lumdis(z, 'conversion'))/u.Mpc)*(3.08567758*10**24)) #conversion to cm
    log_F = row['log_L_agn'] - 2*log_dis - m.log10(4*m.pi)
    return log_F

full['log_F'] = full.apply(lambda row: log_flux(row), axis = 1)


def log_f_l(row):
    log_f_l = row['log_F']
    return log_f_l

full['log_F_l'] = full.apply(lambda row: log_f_l(row), axis = 1)

def dis(row):
    log_d = .5*((m.log10(row['t_cent']) - beta)/alpha + 44 - m.log10(4*m.pi) - row['log_F_l']) #in units of log10 cm
    d = (10**log_d)/(3.08567758*10**24) #converts to linear and cm
    return d

full['dis_mpc'] = full.apply(lambda row: dis(row), axis = 1)
"""
#fix error later. 

#error in just Flux and tau
def dis_Ft_err(row): 
    err_F = 10**(row['log_L_agn_err'])
    err_t = row['t_c_perr']
    F = row['F']
    t = row['t_cent']
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+F
    ddF =  -1*m.sqrt((10**(-1*beta) * t)**(1/alpha))/(4* m.sqrt(m.pi)*F**(3/2))
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+t
    ddt = m.sqrt((10**(-1*beta)* t)**(1/alpha))/(4* m.sqrt(m.pi) * alpha * m.sqrt(F) * t)
    
    dis_err = m.sqrt((ddF*err_F)**2+(ddt*err_t)**2)
    return dis_err
              
full['dis_Ft_err'] = full.apply(lambda row: dis_Ft_err(row), axis = 1)

#error in flux, tau, alpha, and beta
def dis_Ftab_err(row):
    err_a = .0275
    err_b = .024
    F = row['F']
    t = row['t_cent']
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+a
    dda = -1*(m.sqrt((10**(-1*beta) * t)**(1/alpha)) * m.log(10**(-1*beta) * t))/(4 * m.sqrt(m.pi) * alpha**2 * m.sqrt(F))
    #http://www.wolframalpha.com/input/?i=derivative+of+(sqrt((t%2F(10)%5Eb)%5E1%2Fa))%2F(sqrt(4pi*F))+with+respect+to+b
    ddb =  -1*(m.log(10) * m.sqrt((10**(-1*beta) * t)**(1/alpha)))/(4 * m.sqrt(m.pi) * alpha * m.sqrt(F))
    
    dis_err_full = m.sqrt( row['dis_Ft_err']**2 + (dda*err_a)**2 + (ddb*err_b)**2 )
    return dis_err_full 

full['dis_Ftab_err'] = full.apply(lambda row: dis_Ftab_err(row), axis = 1)
"""
def constantmultiple(row):
    point = row['dis_mpc']
    z = row['z']
    real = float((lumdis(z, 's'))/u.Mpc)
    error = real/point
    return m.log10(error)

full['error_order'] = full.apply(lambda row: constantmultiple(row), axis = 1)

print full['error_order'].mean()

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
D_e= 0 #full['dis_Ft_err']
D_ee= 0 #full['dis_Ftab_err']
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

   


legend = ax1.legend(loc="upper left")


for label in legend.get_texts():
    label.set_fontsize('small')


plt.show()


# Calculated luminosity distances are far too high
# Flux is wrong, can't use R should be distance to observer, not R of broad line region