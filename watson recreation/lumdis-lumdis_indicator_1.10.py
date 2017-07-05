import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import LambdaCDM
import astropy.units as u
import math as m 
import pandas as pd 
import  scipy.odr as o



#############################################SETUP###################################################
#to stop the commandline spam
pd.options.mode.chained_assignment = None  # default='warn'

class cosmo(object):
    def __init__(self, H, Om, Ol, Ok):
        self.H = H
        self.Om = Om
        self.Ol = Ol
        self.Ok = Ok
        
cosmology = {"s":cosmo(70, .3089, .6911, 0), 
             "bs0":cosmo(70,.1,.9, 0), 
             "bs1":cosmo(70, 0, 0, 1), 
             "bs2":cosmo(70,.5 , 1, -.5),
             "conversion":cosmo(71, .3,.7,0)  #values gotten from http://adsabs.harvard.edu/abs/2015PASP..127...67B
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



tables = [full, full1]

#################################################ANALISIS#######################################################

#using the 'clean2+ExtCorr' values on page 22 of Bentz et al. 2013
alpha = .549
beta =  1.559   #'K' in paper

"""

something weird is happening

LambdaCDM doens't actually have any code that has curvature. 
I looked at the source files for astropy (the file im looking at called core.py)
(path is anaconda2/lib/python2.7/site-packages/astropy/cosmology)
the FLRW class, that all of the cosmology classes Inherit, doesnt actually initialize self._Ok0
even though it is referenced by a lot of other classes and used in the calculations (that is the value for curvature)
I tried to cahnge core.py, but I havent found a way to make the silgltly reworked lumdis
function to work with it. 

spicificlly, on line 149 of core.py, I changed it to this:
    
    def __init__(self, H0, Om0, Ode0, Ok0, Tcmb0=2.725, Neff=3.04,
                 m_nu=u.Quantity(0.0, u.eV), Ob0=None, name=None):

        self._Ok0 = self.Ok0
        
Im pretty sure this would be enherited correclty, but it doesn't to be. 
If you could figure this shit out, that would be sweet. 

"""

def lumdis(z, key):
    Ho = cosmology[key].H #hubble constant
    mo = cosmology[key].Om #matter
    l = cosmology[key].Ol #dark energyH
    k = cosmology[key].Ok #curvature
    cos = LambdaCDM(H0 = Ho, Om0 = mo, Ode0 = l) #Ok0 = k should work but it doesn't  
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
    if t_cent == 'tau_cent':
        log_d = .5*((m.log10(row[t_cent]*(1+row['redshift'])) - beta)/alpha + 44 - m.log10(4*m.pi) - row['log_F']) #in units of log10 cm
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
    scatter = .036
    d = row['dis_mpc']
    
    ddF = -0.5
    ddr = m.sqrt((0.5/alpha)**2 + scatter)
    
    
    dis_err_log = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    dis_err_full = m.log(10)*d*dis_err_log 
    return dis_err_full 

full['dis_Ftab_err'] = full.apply(lambda row: dis_Ftab_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full1['dis_Ftab_uperr'] = full1.apply(lambda row: dis_Ftab_err(row,0,'tau_cent_uperr'), axis = 1)
full1['dis_Ftab_loerr'] = full1.apply(lambda row: dis_Ftab_err(row,0,'tau_cent_loerr'), axis = 1)

def dis_FTab_averr(row):
    err = (row['dis_Ftab_uperr']+row['dis_Ftab_loerr'])/2
    return err
    
full1['dis_FTab_averr'] = full1.apply(lambda row: dis_FTab_averr(row), axis = 1)


#see https://www.wikiwand.com/en/Log%E2%80%93log_plot
def regression(x, m, b):
    y = x**m + 10**b
    return y 

####linear regression code that isnt done yet####

log_log_model = o.Model(regression)

data = o.RealData([full['z'],full1['redshift']], [full['dis_mpc'],full1['dis_mpc']], sx = 0, sy =  [full['dis_Ftab_err'],full1['dis_FTab_averr']])
"""
odr_output = o.ODR(data, log_log_model, beta0=[0.,1.])

out = odr_output.run()

out.print()
"""

####closeness of fit to model####

#sums the individual values of Chi from the two tables and devides by the number of points
def x_sum(tables, column):
    data = []
    n = []
    for t in tables:
        data.append(t.loc['X_sum',column])
        n.append(len(t[column]))
    
    return sum(data)/sum(n)


#gives one value for x^2 based off of the cosmology you are compairing it with (including Ok in case it gets fixed)
def general_fit(H0, Om, Ol, Ok = 0):
    cozmol =  LambdaCDM(H0 = H0, Om0 = Om, Ode0 = Ol) #again, Ok0 = k should work but it doent
    
    def cosmol_x(row,z,err,cosmology):
        z = row[z] 
        model = float(cosmology.luminosity_distance(z) /u.Mpc)
        x = ((row['dis_mpc'] - model)**2)/row[err]**2
        return x
    
    #makes a column with a name based off of the cosmology you are testing and produces indivitual values for the sum
    #format: 'x^2_H0_OM_OL_OK'
    current_column = 'x^2_%(H0)f_%(OM)f_%(OL)f_%(OK)f' % {'H0': H0,'OM': Om,'OL': Ol,'OK': Ok}
    
    full[current_column] = full.apply(lambda row: cosmol_x(row,'z','dis_Ftab_err', cozmol), axis = 1) 
    full1[current_column] = full1.apply(lambda row: cosmol_x(row,'redshift','dis_FTab_averr', cozmol), axis = 1) 
    
    full.loc['X_sum',current_column] = full[current_column].sum()
    full1.loc['X_sum',current_column] = full1[current_column].sum()
    
    return x_sum([full, full1], current_column)
    


def chi_plot_values(H0 = 70, Om = .3089, Ol = .6911, Ok = 0):
    #add more conditionals later when we figure out how we are going to do a 5D plot. 
    
    #infastructure for 5d:
    #the following is a nested list of lists for each value on the plot, with the form being [X^2, H0, OM,OL,OK]
    
    data = []
    
    if not isinstance(H0, float): 
        for H0 in H0: 
            x2 = general_fit(H0, Om, Ok)
            coordinate = [x2, H0, Om, Ol, Ok]
            data.append(coordinate)
    
    if not isinstance(Om, float): 
        for Om in Om: 
            x2 = general_fit(H0, Om, Ok)
            coordinate = [x2, H0, Om, Ol, Ok]
            data.append(coordinate)
    
    if not isinstance(Ol, float): 
        for Ol in Ol: 
            x2 = general_fit(H0, Om, Ok)
            coordinate = [x2, H0, Om, Ol, Ok]
            data.append(coordinate)
    
    return data

#plotting chi suqred vs different H0
"""
H0 = list(np.arange(40,500,1))

raw_chi_h = chi_plot_values(H0)

chis = []

for i in raw_chi_h:
    chis.append(i[0])

chi_H0 = plt.subplot()
chi_H0.plot(H0, chis)

chi_H0.set_xlabel('Varying H0, everything else constant', fontsize = 14.5)

chi_H0.set_ylabel('resultant Chi Squred', fontsize = 14.5)
"""
#varying Om

Oms = list(np.arange(0,1,.01))

raw_chi_Om = chi_plot_values(70,Oms)

chis_Om = []

for i in raw_chi_Om:
    chis_Om.append(i[0])

chi_Om = plt.subplot()
chi_Om.plot(Oms, chis_Om)

chi_Om.set_xlabel('Varying Om, everything else constant', fontsize = 14.5)

chi_Om.set_ylabel('resultant Chi Squred', fontsize = 14.5)

######################################### MAKING GRAPHS ######################################################

#Everything below is commented out so that I can produce the Chi Plot without having to plotting the other stuff too 

"""
#Labels:

title = 'Luminosity  Distance  vs  Z  For  Different  Cosmologies'
axis_x =  "$Redshift$"
axis_y1 = 'Luminosity Distance $(MPc)$'
axis_y2 = r"$ \frac{\tau}{\sqrt{\lambda  F}}$"

z = np.arange(0,2,.001)
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
        
ax1.plot(z, norm, 'k', label=r"$\Omega_m = .3089$  $\Omega_\lambda = .6911$  $\Omega_k = 0$  (Current Best Cosmolgy)", linewidth = 1)
ax1.plot(z, bs0, 'r:', label=r"$\Omega_m = .1$  $\Omega_\lambda = .9 $  $\Omega_k = 0$", linewidth = 1)
ax1.plot(z, bs1, 'y--', label=r"$\Omega_m = 0$   $\Omega_\lambda = 0$    $\Omega_k = 1$", linewidth = 1)
ax1.plot(z, bs2, 'g-.', label=r"$\Omega_m = .5$  $\Omega_\lambda = 1$  $\Omega_k = -.5$", linewidth = 1)
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
"""