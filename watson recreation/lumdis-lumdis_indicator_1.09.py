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


#regular data

trim = pd.DataFrame()
trim = raw[['name1','log_L_agn','log_L_agn_err','t_cent', 't_c_perr', 'z']]

def percenterror(row):
    return (row['t_c_perr']/row['t_cent'])

trim['pct_err_t'] = trim.apply(lambda row: percenterror(row), axis = 1)

idx =  trim.groupby(['name1'])['pct_err_t'].transform(min) == trim['pct_err_t']
full = trim[idx]

#combination of table1_ascii and mbh_measurments

raw1 = pd.read_csv('mbh_measurements.csv')
raw2 = pd.read_csv('table1_ascii.csv')

raw1.set_index('#RMID')
raw2.set_index('#RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full1 = full1[full1.wLw != -99.0]

#this table has no error in the luminosity, for whatever reason, I dont fucking know. 


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

def log_flux(row,redshift,log_L_agn):
    z = row[redshift]
    log_dis = m.log10(float((lumdis(z, 'conversion'))/u.Mpc)*(3.08567758*10**24)) #conversion to cm
    log_F = row[log_L_agn] - 2*log_dis - m.log10(4*m.pi)
    return log_F

full['log_F'] = full.apply(lambda row: log_flux(row, 'z', 'log_L_agn'), axis = 1)
full1['log_F'] = full1.apply(lambda row: log_flux(row, 'redshift', 'wLw'), axis = 1)

def dis(row,t_cent):
    log_d = .5*((m.log10(row[t_cent]) - beta)/alpha + 44 - m.log10(4*m.pi) - row['log_F']) #in units of log10 cm
    d = (10**log_d)/(3.08567758*10**24) #converts to linear and cm
    return d

full['dis_mpc'] = full.apply(lambda row: dis(row,'t_cent'), axis = 1)
full1['dis_mpc'] = full1.apply(lambda row: dis(row,'tau_cent'), axis = 1)

#error in just Flux and tau
def dis_Ft_err(row, log_L_agn_err,t_c_err): 
    #This is cause there is no L error in one of the data groups
    if log_L_agn_err==0:
        err_F = 0
    else: 
        err_F = row[log_L_agn_err]
    err_r = 1/(m.log(10)*row[t_c_err])
    d = row['dis_mpc']
    
    #https://www.wolframalpha.com/input/?i=derivatibe+of+.5((+r+-b)%2Fa%2B44-log10(4pi)-F)+with+respect+to+F
    ddF = -0.5
    #https://www.wolframalpha.com/input/?i=derivatibe+of+.5((+r+-b)%2Fa%2B44-log10(4pi)-F)+with+respect+to+r
    ddr = 0.5/alpha
    
    dis_err_log = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    #https://www.wolframalpha.com/input/?i=derivative+10%5Ex
    dis_err = m.log(10)*d*dis_err_log #because the error in d is log but d itself is linear, 10^log10(d) just produces d 
    return dis_err
              
full['dis_Ft_err'] = full.apply(lambda row: dis_Ft_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full1['dis_Ft_uperr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_uperr'), axis = 1)
full1['dis_Ft_loerr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_loerr'), axis = 1)

#error in flux, tau, alpha, and beta NOW USING SCATTER!
def dis_Ftab_err(row, log_L_agn_err,t_c_err):
    if log_L_agn_err==0:
        err_F = 0
    else: 
        err_F = row[log_L_agn_err]
    err_r = 1/(m.log(10)*row[t_c_err])
    #like aplha and beta, it comes from clean2 + ExtCorr
    scatter = .024
    d = row['dis_mpc']
    
    ddF = -0.5
    ddr = m.sqrt((0.5/alpha)**2 + scatter)
    
    
    dis_err_log = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    dis_err_full = m.log(10)*d*dis_err_log 
    return dis_err_full 

full['dis_Ftab_err'] = full.apply(lambda row: dis_Ftab_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full1['dis_Ftab_uperr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_uperr'), axis = 1)
full1['dis_Ftab_loerr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_loerr'), axis = 1)

def constantmultiple(row):
    point = row['dis_mpc']
    z = row['z']
    real = float((lumdis(z, 's'))/u.Mpc)
    error = real/point
    return m.log10(error)

full['error_order'] = full.apply(lambda row: constantmultiple(row), axis = 1)

print full
print full1

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

#from data 1
Z = full['z']
D = full['dis_mpc']
D_e= full['dis_Ft_err']
D_ee= full['dis_Ftab_err'] 

#from data 2
Z1 = full1['redshift']
D1 = full1['dis_mpc']
D_e1= [full1['dis_Ft_loerr'],full1['dis_Ft_uperr']]
D_ee1= [full1['dis_Ftab_loerr'],full1['dis_Ftab_uperr']] 
        
ax1.loglog(z, norm, 'k', label=r"$\Omega_m = .3089$  $\Omega_\lambda = .6911$  $\Omega_k = 0$  (Current Best Cosmolgy)", linewidth = 1)
ax1.loglog(z, bs0, 'r:', label=r"$\Omega_m = .1$  $\Omega_\lambda = .9 $  $\Omega_k = 0$", linewidth = 1)
ax1.loglog(z, bs1, 'y--', label=r"$\Omega_m = 0$   $\Omega_\lambda = 0$    $\Omega_k = 1$", linewidth = 1)
ax1.loglog(z, bs2, 'g-.', label=r"$\Omega_m = .5$  $\Omega_\lambda = 1$  $\Omega_k = -.5$", linewidth = 1)
ax1.set_xlabel(axis_x, fontsize = 14.5)
plt.title(title, fontsize = 15)
ax1.set_ylabel(axis_y1, fontsize = 14.5)
#set 1
ax1.errorbar(Z, D,xerr=0, yerr = D_e, fmt = 'k.', ecolor='r'  , elinewidth= 1.4, label=r'Redshift Independant Distance $without$ $\alpha$ and $\beta$ error',zorder=10)
ax1.errorbar(Z, D,xerr=0, yerr = D_ee, fmt = 'k.', ecolor='#FF9900'  , elinewidth= 1.3, label=r'Redshift Independant Distance $with$ $\alpha$ and $\beta$ error',zorder=9)
#set 2
ax1.errorbar(Z1, D1,xerr=0, yerr = D_e1, fmt = 'b.', ecolor='r'  , elinewidth= 1.4, label=r'Redshift Independant Distance $without$ $\alpha$ and $\beta$ error',zorder=10)
ax1.errorbar(Z1, D1,xerr=0, yerr = D_ee1, fmt = 'b.', ecolor='#FF9900'  , elinewidth= 1.3, label=r'Redshift Independant Distance $with$ $\alpha$ and $\beta$ error',zorder=9)

#ax1.axis([10**-3,4,4,1.1*10**5])


legend = ax1.legend(loc="upper left")


for label in legend.get_texts():
    label.set_fontsize('small')

plt.show()