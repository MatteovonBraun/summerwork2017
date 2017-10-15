'''
=====================================================================================================

Montey Carlo Simulation of the Luminsity - Radius Relationship
  
Version: 0.2

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
 
- Get rid of duplcates 
- Add in other data (make sure the error and everything else is propperly logged. )
    . there is positive and negative error in some, so just look for what I did in 
    . the other file. 
- Add the total implementation of what is in Yasaman's Email
- Make a proper distribution of luminosity of QSO's to go in the Yasaman code. 
- Make everything a function so that it can be used by the main file. 

"""


############# INITIALIZATION ##############
raw1 = pd.read_csv('mbh_measurements.csv')
raw2 = pd.read_csv('table1_ascii.csv')

raw1.set_index('#RMID')
raw2.set_index('#RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full1 = full1[full1.wLw != -99.0]

############# COMPUTATION ###########################

def t_log(row, column):
    log = m.log10(row[column])
    return log 
    
    
full1['log_t_cent'] = full1.apply(lambda row: t_log(row,'tau_cent'), axis = 1)

def t_err_log(row, val_column, err_column):
    err = row[val_column]/(row[err_column]*m.log(10))
    return err


full1['log_t_c_perr'] = full1.apply(lambda row: t_err_log(row,'tau_cent','tau_c_perr'), axis = 1)

#to deal with that '10^44' in equation (2) in Bentz et. al. 
def L_44(row):
    return row['log_L_agn'] - 44
    
full1['log_L_agn_44'] = full1.apply(lambda row: L_44(row), axis = 1)

#LINE OF BEST FIT#

L = full1['log_L_agn_44']
t = full1['log_t_cent']
L_err = full1['log_L_agn_err']
t_err = full1['log_t_c_perr'] 

#remember, everything is in log, so when we are converting back we will have 
#to convert the linear functio too. 
#again:, 
#https://stackoverflow.com/questions/22670057/linear-fitting-in-python-with-uncertainty-in-both-x-and-y-coordinates


#WE HAVE TO CHANGE THIS TO MULTIPLE LINEAR REGRESSION


def linear(p, x):
    m, c = p
    return m*x + c

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

x0 = list(np.arange(min(full1['log_L_agn_44']),max(full1['log_L_agn_44']),.01))
y0 = list(map(lambda x: x*m+c, x0))


a = .549
k = 1.559 
x1 = x0
y1 = list(map(lambda x: x*a+k, x1))

pl.errorbar(L, t,xerr=L_err, yerr = t_err, fmt = 'k.', ecolor='r'  , elinewidth= .5, label=r'Log L vs Log $\tau$',zorder=10)
pl.plot(x0,y0,'b:', label='Best fit using ODR of this data', linewidth=1.2)
pl.plot(x1,y1,'g--',label='Best fit using Bentz et. al. numbers',linewidth=.6)

legend = pl.legend(loc="upper left")
pl.set_xlabel('Log Luminosity', fontsize = 14.5)
pl.set_ylabel('Log Tau', fontsize = 14.5)


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
