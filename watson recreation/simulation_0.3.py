'''
=====================================================================================================

Montey Carlo Simulation of the Luminsity - Radius Relationship
  
Version: 0.1

Authors: McDonnell, Matteo; Moore, Marc

=====================================================================================================
'''
import numpy as np
import math as m
import matplotlib.pyplot as plt
import pandas as pd 
import scipy.odr as o

"""
TO DO:
   
- Add the total implementation of what is in Yasaman's Email

- Make a proper distribution of luminosity of QSO's to go in the Yasaman code. 

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

raw1['name1'] = raw1['RMID']
raw2['name1'] = raw2['RMID']

raw1.set_index('RMID')
raw2.set_index('RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full1 = full1[full1.wLw != -99.0]


def aveerr(row):
    return (row['tau_cent_uperr']+row['tau_cent_loerr'])/2

full1['t_c_err'] = full1.apply(lambda row: aveerr(row), axis = 1)

full1 = full1.rename(columns = {'redshift':"z",'wLw':'log_L_agn'})

def fix_tau(row):
    return row['tau_cent']*(row['z']+1)
    
full1['t_cent'] = full1.apply(lambda row: fix_tau(row), axis = 1)    


#This is so we can add the two dataframes cleanly and not have a devide by 0 error in ODR
full1['log_L_agn_err'] = 0.000000000000001

#data 1 + 2

prep0 = full0[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]
prep1 = full1[['log_L_agn','log_L_agn_err','t_cent','t_c_err','z']]

full2 = pd.concat([prep0,prep1])


############# COMPUTATION ###########################


def t_log(row, column):
    log = m.log10(row[column])
    return log 
    
    
full0['log_t_cent'] = full0.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full1['log_t_cent'] = full1.apply(lambda row: t_log(row,'t_cent'), axis = 1)
full2['log_t_cent'] = full2.apply(lambda row: t_log(row,'t_cent'), axis = 1)

def t_err_log(row, val_column, err_column):
    err = row[err_column]/(row[val_column]*m.log(10))
    return err

full0['log_t_c_perr'] = full0.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full1['log_t_c_perr'] = full1.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)
full2['log_t_c_perr'] = full2.apply(lambda row: t_err_log(row,'t_cent','t_c_err'), axis = 1)

#to deal with that '10^44' in equation (2) in Bentz et. al. 
def L_44(row):
    return row['log_L_agn'] - 44
    
full0['log_L_agn_44'] = full0.apply(lambda row: L_44(row), axis = 1)
full1['log_L_agn_44'] = full1.apply(lambda row: L_44(row), axis = 1)
full2['log_L_agn_44'] = full2.apply(lambda row: L_44(row), axis = 1)


#LINE OF BEST FIT#

L = full2['log_L_agn_44']
t = full2['log_t_cent']
L_err = full2['log_L_agn_err']
t_err = full2['log_t_c_perr'] 

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

############## PLOT-ATION ################


fig, pl = plt.subplots()

#using the output of ODR:

x0 = list(np.arange(min(full0['log_L_agn_44']),max(full0['log_L_agn_44']),.01))
y0 = list(map(lambda x: x*m+c, x0))

a = .549
k = 1.559 
x1 = x0
y1 = list(map(lambda x: x*a+k, x1))

#datapoints
L0 = full0['log_L_agn_44']
t0 = full0['log_t_cent']
L_err0 = full0['log_L_agn_err']
t_err0 = full0['log_t_c_perr'] 

L1 = full1['log_L_agn_44']
t1 = full1['log_t_cent']
L_err1 = full1['log_L_agn_err']
t_err1 = full1['log_t_c_perr'] 


pl.errorbar(L0, t0,xerr=L_err0, yerr = t_err0, fmt = 'k.', ecolor='r'  , elinewidth= .5, label=r'Log L vs Log $\tau$ Misty Data',zorder=10)
pl.errorbar(L1, t1,xerr=L_err1, yerr = t_err1, fmt = 'b.', ecolor='r'  , elinewidth= .5, label=r'Log L vs Log $\tau$ New data (source?)',zorder=10)
pl.plot(x0,y0,'b:', label='Best fit using ODR of this data', linewidth=1.2)
pl.plot(x1,y1,'g--',label='Best fit using Bentz et. al. numbers',linewidth=.6)

legend = pl.legend(loc="upper left")
pl.set_xlabel('Log Luminosity', fontsize = 14.5)
pl.set_ylabel('Log Tau', fontsize = 14.5)

#pl.axis([-2.6,2,0,2.5])

plt.show()


"""

THIS IS IDL CODE THAT JON SENT 

THIS IS EVENTUALLY WHAT WE WILL BE DOING FOR THE DISTROBUTION OF QSOS



  nsim = 1000   ;number of simulated points
  minlogL = 42   ;range of log(L)
  maxlogL = 47
  logL = findgen(nsim) / (nsim-1) * (maxlogL-minlogL) + minlogL

  alpha = -3.4   ;constants for luminosity function
  beta = -1.6
  Mbreak = -21
  logLbreak = -0.4 * (Mbreak - 88.7)

  L = 10^(logL/logLbreak)   ;easier parameterization
  P = 1 / (L^alpha + L^beta)  ;luminosity function

  bad = lindgen(nsim)
  Parr = fltarr(nsim)
  logLdist = fltarr(nsim)
  yarr = fltarr(nsim)
  nbad = nsim
  niter = 0

  repeat begin
     logLdist[bad] = randomu(seed,nbad) * (maxlogL-minlogL) + minlogL
     yarr[bad] = randomu(seed,nbad)

     L = 10^(logLdist / logLbreak)
     parr[bad] = 1 / (L[bad]^alpha + L[bad]^beta)

     bad = where(yarr gt parr,nbad)
     niter++
;;      plothist,Ldist,/ylog,bin=0.1,xrange=[minlogL,maxlogL]  ;check distribution
;;      djs_oplot,logL,P*0.05*nsim/max(P),color='red'  ;target distribution in red
;;      print,niter,nbad
  endrep until nbad eq 0

  return,logLdist
"""