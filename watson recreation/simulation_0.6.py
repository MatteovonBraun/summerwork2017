'''
=====================================================================================================

Montey Carlo Simulation of the Luminsity - Radius Relationship
  
Version: 0.6

Authors: McDonnell, Matteo; Moore, Marc

=====================================================================================================
'''
import numpy as np
import math as m
import matplotlib.pyplot as plt
import pandas as pd 
import scipy.odr as o
import scipy.stats as st
from scipy.integrate import quad
from decimal import *


"""
TO DO:
   
- See if you just simply add the scatter and the random error in the simulation (see simgen)

- Use a different way of dealing with the numbers in the real distrobution so
  that it actuall works? (see the usage of 'decimal')

- Finish the rest of Monte Carlo 

- Make everything a function so that it can be used by the main file. 

"""

#to stop the commandline spam
pd.options.mode.chained_assignment = None  # default='warn'

#to set the precsision of the decimal package
getcontext().prec = 35
getcontext().Emax = 9999999999999999999999

############# INITIALIZATION ##############

#data 0
raw0 = pd.read_csv('CSVResults z.csv')
raw0.sort_values(by='name1')


trim = raw0[['name1','log_L_agn','log_L_agn_err','t_cent', 't_c_perr', 'z']]

def percenterror(row):
    return (row['t_c_perr']/row['t_cent'])

trim['pct_err_t'] = trim.apply(lambda row: percenterror(row), axis = 1)

idx =  trim.groupby(['name1'])['pct_err_t'].transform(min) == trim['pct_err_t']
full0 = trim[idx]

full0 = full0.rename(columns={'t_c_perr':'t_c_err'})

full0.reset_index

full0.drop('pct_err_t', axis=1, inplace=True)

#data 1

raw1 = pd.read_csv('mbh_measurements.csv')
raw2 = pd.read_csv('table1_ascii.csv')

raw1['name1'] = raw1['RMID']
raw2['name1'] = raw2['RMID']

raw1.set_index('RMID')
raw2.set_index('RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full1 = full1[full1.wLw != -99.0]


def aveerr(row, up, down):
    return (row[up]+row[down])/2

full1['t_c_err'] = full1.apply(lambda row: aveerr(row,'tau_cent_uperr','tau_cent_loerr'), axis = 1)

full1 = full1.rename(columns = {'redshift':"z",'wLw':'log_L_agn'})

def fix_tau(row):
    return row['tau_cent']*(row['z']+1)
    
full1['t_cent'] = full1.apply(lambda row: fix_tau(row), axis = 1)    


#This is so we can add the two dataframes cleanly and not have a devide by 0 error in ODR
full1['log_L_agn_err'] = 0.2

#data 3

raw3 = pd.read_csv('du16_with_z.csv')

#raw3['name1'] = raw3['Object']

raw3['t_c_err'] = raw3.apply(lambda row: aveerr(row,'t_c_perr','t_c_nerr'), axis = 1)
full2 = raw3


#data 1 + 2 +3

prep0 = full0[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]
prep1 = full1[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]
prep2 = full2[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]


full3 = pd.concat([prep0,prep1,prep2])


####################### COMPUTATION ###########################


def t_log(row, column):
    log = m.log10(row[column])
    return log 
    
    
full0['log_t_cent'] = full0.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full1['log_t_cent'] = full1.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full2['log_t_cent'] = full2.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full3['log_t_cent'] = full3.apply(lambda row: t_log(row,'t_cent'), axis = 1)

def t_err_log(row, val_column, err_column):
    err = row[err_column]/(row[val_column]*m.log(10))
    return err

full0['log_t_c_err'] = full0.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full1['log_t_c_err'] = full1.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full2['log_t_c_err'] = full2.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full3['log_t_c_err'] = full3.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)

#to deal with that '10^44' in equation (2) in Bentz et. al. 
def L_44(row):
    return row['log_L_agn'] - 44
    
full0['log_L_agn_44'] = full0.apply(lambda row: L_44(row), axis = 1)
full1['log_L_agn_44'] = full1.apply(lambda row: L_44(row), axis = 1)
full2['log_L_agn_44'] = full2.apply(lambda row: L_44(row), axis = 1)
full3['log_L_agn_44'] = full3.apply(lambda row: L_44(row), axis = 1)


#LINE OF BEST FIT#

L = full3['log_L_agn_44']
t = full3['log_t_cent']
L_err = full3['log_L_agn_err']
t_err = full3['log_t_c_err'] 

#remember, everything is in log, so when we are converting back we will have 
#to convert the linear functio too. 
#again:, 
#https://stackoverflow.com/questions/22670057/linear-fitting-in-python-with-uncertainty-in-both-x-and-y-coordinates


def linear(p, x):
    
    return p[0]*x + p[1]

linear_model = o.Model(linear)

data = o.RealData(L, t, sx = L_err, sy = t_err)

odr = o.ODR(data, linear_model,beta0 = [0.,1.])

output = odr.run() 

m = output.beta[0]
c = output.beta[1]


output.pprint()

print output.beta

################### SIMULATION ##################

#constants of simulation:

"""
Guide for the constants for simulation: 

The first item is the number of objects to be generated 

The first dictionary is for the error to be added to luminosity, 
second for the errors and pertibation in tau. 

intr_er = 'intrinsic scatter in the error'
scat_er_max = 'the maximum ammount of error to be added to the scatter (will only be positive)'
scat_er_sig = 'standard deviation (sigma) of the normal distribution of the added error'

pert_max = 'maximum pertibation of tau up or down (will be both positive and negative'
pert_sig = 'standard deviation of tau pertiba'

The third dictionary  is for the distribution generation, and the names are self explanitory

The fourth dictionary defines the minimum and maximum of the distrobution 
"""

def tau_lum(lum,a,b): #10^44
    return b + a*lum 

scat_l = 1.
scat_t = tau_lum(scat_l,m,c)

sim_c = [10000,
         {'intr_er':scat_l,
          'scat_er_max':1.,
          'scat_er_sig':1.},
         {'intr_er':scat_t,
          'scat_er_max':1.,
          'scat_er_sig':1.,
          'pert_max':1.,
          'pert_sig':1.},
         {'breakL':43.88,
          'alpha':3.4,
          'beta':1.6},
         {'min':42,
          'max':47}]

"""
This doesn't work yet, just using a uniform distrobution 

def dist():

    def pdf(l):
    
        L = (10**Decimal(l/sim_c[3]['breakL']))
        D = 1/(L**Decimal(sim_c[3]['beta']) + (L**Decimal(sim_c[3]['alpha'])))
        return D
    
    correction = quad(pdf, sim_c[4]['min'], sim_c[4]['max'])[0]

    class my_pdf(st.rv_continuous):
        
        def _pdf(self,l_L): 
            #"l_L" in this is always log L     
            
            return float(pdf(l_L))/correction 
            
    dist_Log_L = my_pdf(momtype = 0, a = sim_c[4]['min'] ,b=sim_c[4]['max'],name='l_L_dist')

    distro = dist_Log_L.rvs(size = sim_c[0])

    return distro



def hist_plot(dis,bn):
        
    fig, hst = plt.subplots()
    
    hst.hist(dis, bins = bn)
    
    plt.show()

hist_plot(dist(),20)
"""

#this isn't done, it took me a few goes to format this correctly. Its mostly there. 
def simgen():

    def norm_val(mini, maxi, sigma):

        lower, upper = mini, maxi
        mu, sigma = 0, sigma 
        X = st.truncnorm( (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma )
        i = X.rvs(1)
        return i 

    def rawdist(): 
        #makeing a uniform distrobution (use the above funciton when it acutually works)
        a = sim_c[4]['min']
        b = sim_c[4]['max']
        n = sim_c[0]
        d = np.random.uniform(a,b,n)
        return d 

    def err(n): #1 for luminosity, 2 for tau. 
        scat = sim_c[n]['intr_er']
        rand = norm_val(0,sim_c[n]['scat_er_max'],sim_c[n]['scat_er_sig'])
        tot_err = scat+rand

    def frame():
        draft = pd.DataFrame()
        draft['log_L_agn'] = rawdist()
        






############################ PLOT-ATION ##########################

def L_tau_plot():

    fig, pl = plt.subplots()

    #using the output of ODR:

    x0 = list(np.arange(-3,3,.01))
    y0 = list(map(lambda x: x*m+c, x0))

    a = .549
    k = 1.559 
    x1 = x0
    y1 = list(map(lambda x: x*a+k, x1))


    #datapoints
    L0 = full0['log_L_agn_44']
    t0 = full0['log_t_cent']
    L_err0 = full0['log_L_agn_err']
    t_err0 = full0['log_t_c_err'] 

    L1 = full1['log_L_agn_44']
    t1 = full1['log_t_cent']
    L_err1 = full1['log_L_agn_err']
    t_err1 = full1['log_t_c_err'] 

    L2 = full2['log_L_agn_44']
    t2 = full2['log_t_cent']
    L_err2 = full2['log_L_agn_err']
    t_err2 = full2['log_t_c_err'] 

    #plot labels
    Label0 = r'Log L vs Log $\tau$ Misty Data'
    Label1 = r'Log L vs Log $\tau$ New data (source?)'
    Label2 = r'Log L vs Log $\tau$ Du data'
    Label3 = 'Best fit using ODR of all data'
    Label4 = 'Best fit using Bentz et. al. numbers'

    pl.errorbar(L0, t0,xerr=L_err0, yerr = t_err0, fmt = 'k.', ecolor='k'  , elinewidth= .5, label=Label0,zorder=10)
    pl.errorbar(L1, t1,xerr=L_err1, yerr = t_err1, fmt = 'b.', ecolor='b'  , elinewidth= .5, label=Label1,zorder=10)
    pl.errorbar(L2, t2,xerr=L_err2, yerr = t_err2, fmt = 'g.', ecolor='g'  , elinewidth= .5, label=Label2,zorder=10)
    pl.plot(x0,y0,'b:', label=Label3, linewidth=1.2)
    pl.plot(x1,y1,'g--',label=Label4, linewidth= .6)

    legend = pl.legend(loc="upper left")
    pl.set_xlabel('Log Luminosity', fontsize = 14.5)
    pl.set_ylabel('Log Tau', fontsize = 14.5)

    pl.axis([-3,3,0,3])

    plt.show()

L_tau_plot()