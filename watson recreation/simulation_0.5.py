'''
=====================================================================================================

Montey Carlo Simulation of the Luminsity - Radius Relationship
  
Version: 0.5

Authors: McDonnell, Matteo; Moore, Marc

=====================================================================================================
'''
import numpy as np
import math as m
import matplotlib.pyplot as plt
import pandas as pd 
import scipy.odr as o
import scipy.stats as st

"""
TO DO:
   
- Add the total implementa
ion of what is in Yasaman's Email

- Make a proper distribution of luminosity of QSO's to go in the Yasaman code. 
    . You can use rv_continous hopefully? IDK how it works, the scipy website
    . isn't working at the time of writing this. 

- Make everything a function so that it can be used by the main file. 

"""

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
raw3 = pd.read_csv('du16.csv')

raw1['name1'] = raw1['RMID']
raw2['name1'] = raw2['RMID']
raw3['name1'] = raw3['NAME']

raw1.set_index('RMID')
raw2.set_index('RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]
trim3 = raw3[['tau(Hb)','errhi','errlo']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full2 = pd.DataFrame(full2.join(trim3, how='inner'))
full2 = full2[full2.wLw != -99.0]


def aveerr(row):
    return (row['tau_cent_uperr']+row['tau_cent_loerr'])/2

full2['t_c_err'] = full2.apply(lambda row: aveerr(row), axis = 1)

full2 = full2.rename(columns = {'redshift':"z",'wLw':'log_L_agn'})

def fix_tau(row):
    return row['tau_cent']*(row['z']+1)
    
full2['t_cent'] = full2.apply(lambda row: fix_tau(row), axis = 1)    


#This is so we can add the two dataframes cleanly and not have a devide by 0 error in ODR
full2['log_L_agn_err'] = 0.2

#data 1 + 2

prep0 = full0[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]
prep1 = full2[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]

full2 = pd.concat([prep0,prep1])


####################### COMPUTATION ###########################


def t_log(row, column):
    log = m.log10(row[column])
    return log 
    
    
full0['log_t_cent'] = full0.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full2['log_t_cent'] = full2.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full2['log_t_cent'] = full2.apply(lambda row: t_log(row,'t_cent'), axis = 1)

def t_err_log(row, val_column, err_column):
    err = row[err_column]/(row[val_column]*m.log(10))
    return err

full0['log_t_c_err'] = full0.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full2['log_t_c_err'] = full2.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full2['log_t_c_err'] = full2.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)

#to deal with that '10^44' in equation (2) in Bentz et. al. 
def L_44(row):
    return row['log_L_agn'] - 44
    
full0['log_L_agn_44'] = full0.apply(lambda row: L_44(row), axis = 1)
full2['log_L_agn_44'] = full2.apply(lambda row: L_44(row), axis = 1)
full2['log_L_agn_44'] = full2.apply(lambda row: L_44(row), axis = 1)


#LINE OF BEST FIT#

L = full2['log_L_agn_44']
t = full2['log_t_cent']
L_err = full2['log_L_agn_err']
t_err = full2['log_t_c_err'] 

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
n_obj = 10000

"""
Guide for the constants for simulation: 

The first list is for the error to be added to luminosity, 
second for the errors and pertibation in tau. 

intr_er = 'intrinsic scatter in the error'
scat_er_max = 'the maximum ammount of error to be added to the scatter (will only be positive)'
scat_er_sig = 'standard deviation (sigma) of the normal distribution of the added error'

pert_max = 'maximum pertibation of tau up or down (will be both positive and negative'
pert_sig = 'standard deviation of tau pertiba'
"""

def tau_lum(lum,a,b): #10^44
    return b + a*lum 

scat_l = 1.
scat_t = tau_lum(scat_l,m,c)

sim_c = [{'intr_er':scat_l,
          'scat_er_max':1.,
          'scat_er_sig':1.},
         {'intr_er':scat_t,
          'scat_er_max':1.,
          'scat_er_sig':1.,
          'pert_max':1.,
          'pert_sig':1.}]




#Distribution of Luminositys (change to real distrobution later)
#These are values that need to be multiplied by 10^44 at the end

dist = st.uniform()



"""
I am pretty sure that this will work with the production of a
proper distrobution of luminisity of QSO's. Im not sure about
the constants, but I'm pretty sure the function is right.

logLbreak = 43.88 #change these later
alpha = 3.4
beta = 1.6


class my_pdf(st.rv_continuous):
    
    def _pdf(self,l_L): 
        #"l_L" in this is always log L        
        L = 10**(l_L/logLbreak)
        D = 1/(L**alpha + L**beta)
        return D
        
dist_Log_L = my_pdf(momtype = 0, a = 0,name='l_L_dist')



#Dont do this, there is an easyer way (just do dist = dist_Log_L.rvs(size = [size of dataset])
def dist_gen(itt):
    dist = []   
    
    print "Generating",itt,'datapoints...'
    print
    
    for i in range(itt):
        dist.append(dist_Log_L.rvs())
        perc = float(i)*100/float(itt)                  
        if (perc%5 == 0):
            
            print  str(round(perc,1)).zfill(3),'% done'

    print '100  % done'       

    return dist

def hist_plot(bn,itt):
    distro = dist_gen(itt)
    
    fig, hst = plt.subplots()
    
    hst.hist(distro, bins = bn)
    
    plt.show()
    
hist_plot(20, 3000)

"""

############################ PLOT-ATION ##########################


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

L1 = full2['log_L_agn_44']
t1 = full2['log_t_cent']
L_err1 = full2['log_L_agn_err']
t_err1 = full2['log_t_c_err'] 

#plot labels
Label0 = r'Log L vs Log $\tau$ Misty Data'
Label1 = r'Log L vs Log $\tau$ New data (source?)'
Label2 = 'Best fit using ODR of this data'
Label3 = 'Best fit using Bentz et. al. numbers'

pl.errorbar(L0, t0,xerr=L_err0, yerr = t_err0, fmt = 'k.', ecolor='r'  , elinewidth= .5, label=Label0,zorder=10)
pl.errorbar(L1, t1,xerr=L_err1, yerr = t_err1, fmt = 'b.', ecolor='r'  , elinewidth= .5, label=Label1,zorder=10)
pl.plot(x0,y0,'b:', label=Label2, linewidth=1.2)
pl.plot(x1,y1,'g--',label=Label3, linewidth= .6)

legend = pl.legend(loc="upper left")
pl.set_xlabel('Log Luminosity', fontsize = 14.5)
pl.set_ylabel('Log Tau', fontsize = 14.5)

pl.axis([-3,3,0,3])

plt.show()
