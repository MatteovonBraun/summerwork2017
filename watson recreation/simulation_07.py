'''
=====================================================================================================

Montey Carlo Simulation of the Luminsity - Radius Relationship
  
Version: 0.6

Authors: McDonnell, Matteo; Moore, Marc

=====================================================================================================
'''
import numpy as np
import math as m
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd 
import scipy.odr as o
import scipy.stats as st
from scipy.integrate import quad
from decimal import *
import datetime as dt


"""
TO DO:

- Make everything much more general

- Make everything a function so that it can be used by the main file. 

"""

#to stop the commandline spam
pd.options.mode.chained_assignment = None  # default='warn'

#to set the precsision of the decimal package
getcontext().prec = 35
getcontext().Emax = 9999999999999999999999

startTime = dt.datetime.now()

def get_time():
    return (dt.datetime.now() - startTime)


############# INITIALIZATION ##############

"""
This is so we dont have to spam .apply functions everywhere. 

You pass the funciton name to apply,

a list of lists columns, where the first item in the sublist is the new 
column name and the second item is a list of parameters, 

and a list of lists dataframes to apply this function to, where the 
first item in the sublist is the frame to make the column in, and 
the second is the the frame on which the function is applied.  
"""

def func_app(function,columns,frames):

    for column in columns: 

        if len(column) == 1:

            for frame in frames: 

                newframe = frame[0]
                oldframe = frame[1]
                newframe[column[0]] = oldframe.apply(lambda row: function(row), axis = 1)

        else: 

            for frame in frames:

                newframe = frame[0]
                oldframe = frame[1]
                newframe[column[0]] = oldframe.apply(lambda row: function(row,column[1]), axis = 1)


#this is to allow the list of inputs to be put into a dicitonary so you can more easilly see the name of a constant. 
def input_unpack(item_list,inputs):

    constants = {}
    
    for n in range(len(inputs)):

        constants.update( { item_list[n] : inputs[n] } )

    return constants


def aveerr(row, inputs):

    c_ = input_unpack(['lo_err','up_err'],inputs)

    return (row[c_['lo_err']]+row[c_['up_err']])/2



#data 1
full1 = pd.read_csv('Bentz_formatted_clean.csv')
full1.sort_values(by='Object')


#data 2

full2 = pd.read_csv('SDSS_new_t.csv')

#data 3

full3 = pd.read_csv('du16_with_z.csv')

#data 1 + 2 + 3
full0 = pd.concat([full1,full2,full3])
full0 = full0.reset_index(drop=True)

frames_ = [[full0,full0],[full1,full1],[full2,full2],[full3,full3]]

###################################


func_app(aveerr,[['t_c_aerr',['t_c_nerr','t_c_perr']]], frames_)


####################### COMPUTATION ###########################

def t_log(row, inputs):

    c_ = input_unpack(['up_err'],inputs)

    log = m.log10(row[c_['up_err']])
    return log 
    

def t_err_log(row, inputs):
    
    c_ = input_unpack(['val','err'],inputs)

    err = row[c_['err']]/(row[c_['val']]*m.log(10))
    return err


#to deal with that '10^44' in equation (2) in Bentz et. al. 
def L_44(row):
    return row['log_L_agn'] - 44
    

func_app(t_log,[['log_t_cent',['t_cent']]],frames_)

func_app(t_err_log,[['log_t_c_aerr',['t_cent','t_c_aerr']],
                    ['log_t_c_nerr',['t_cent','t_c_nerr']],
                    ['log_t_c_perr',['t_cent','t_c_perr']]],frames_)

func_app(L_44,[['log_L_agn_44']],frames_)




#show the highest value of postive log tau error
#print full2.ix[full2['log_t_c_perr'].idxmax()]


#LINE OF BEST FIT#


def linear_regression(): 

    frame = full0

    L = frame['log_L_agn_44']
    t = frame['log_t_cent']
    L_err = frame['log_L_agn_err']
    t_err = frame['log_t_c_aerr'] 

    #remember, everything is in log, so when we are converting back we will have 
    #to convert the linear function too. 
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
    print 
    print output.beta
    print 

    return m,c

m,c = linear_regression()


################### SIMULATION ##################

#constants of simulation:


def frame_mode(mode):
    #0 for all data, 1 for Bentz, 2 for SDSSRM, 3 for Du
    
    frame_ = frames_[mode]
    return frame_


#constants for the distribution (dicitonary 3) are from the fourth row of table 2 in Boyle et al. 2000
sim_mode = 3

sim_n = 100

set_ = frame_mode(sim_mode)
set_0 = set_[0]



sim_c = {'breakL0':44.62,
        'alpha':3.37,
        'beta':1.16,
        'k1':1.241,
        'k2':-.249,
        'l_min':set_0['log_L_agn'].min(),
        'l_max':set_0['log_L_agn'].max(),
        'z_min':set_0['z'].min(),
        'z_max':set_0['z'].max(),
        't_min':set_0['log_t_cent'].min(),
        't_max':set_0['log_t_cent'].max()}


convert_scatter = .018


#this returns both a boolian if a percentage is at a certain value, and the percentage itself (used in console printout)
def bool_perc(n,tot):
    perc_interval = 5
    bool_ = ((float(n)/float(tot)*100)%perc_interval == 0) and (n != 0)
    perc_ = str(int(float(n)/float(tot)*100)).zfill(2)
    return bool_,perc_

def tau_lum(lum,a,b): #10^44
    return float(b + a*lum)

#this defines the luminosity distrobution. If it is used with a 
def dist(row):

    k1 = sim_c['k1']
    k2 = sim_c['k2']


    #this is for the printout in terminal
    if bool_perc(row.name,sim_n)[0]: 
        print bool_perc(row.name,sim_n)[1]+'%  ... ',get_time(), 'has passed.'
    
    z = row['z']
    breakL = sim_c['breakL0']*10**(k1*z + k2*z**2)
    size_ = 1  

    def pdf(l):

        L = Decimal(10**(l-breakL))
        D = 1/(L**Decimal(sim_c['beta']-1) + (L**Decimal(sim_c['alpha']-1)))
        return float(D)
    
    correction = quad(pdf, sim_c['l_min'], 100)[0]
   
    class my_pdf(st.rv_continuous):
        
        def _pdf(self,l_L): 
            #"l_L" in this is always log L     
            return float(pdf(l_L))/correction 
            
    dist_Log_L = my_pdf(momtype = 0, a = sim_c['l_min'],name='l_L_dist')

    distro = dist_Log_L.rvs(size = size_)
  
    return distro



#for testing purposes 
def hist_plot(dis,bn):
        
    fig, hst = plt.subplots()
    
    hst.hist(dis,bins = bn)
    #plt.yscale('log')
    plt.show()

#hist_plot(dist(),8)


def simgen(bounds):

    frame_ = frame_mode(bounds)


    def redshift_dist(): #returns a flat distribution of size sim_n


        z = st.uniform(sim_c['z_min'], sim_c['z_max'])
        z_dist = z.rvs(size = sim_n)
        return z_dist 

    def norm_val(mini, maxi, sigma): #returns a single value that is generated from a truncated normal distrobution

        lower, upper = mini, maxi
        mu, sigma = 0, sigma 
        X = st.truncnorm( (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma )
        i = X.rvs(1)
        return float(i) 

    def l_err(row): #1 for luminosity, 2 for tau. Returns a value of error that adds together 
        
        column = 'log_L_agn_err'
        place = int(np.random.randint(len(set_0.index), size=1))

        sig_ = set_0[column][place]
        max_ = 5 * sig_
        rand = float(norm_val(0,max_,sig_))
        tot_err = rand
        return tot_err

    def t_err(row):


        column = 'log_L_agn_err'
        place = int(np.random.randint(len(set_0.index), size=1))
        scat = convert_scatter #intrinsic scatter in the conversion 

        sig_ = set_0[column][place]
        max_ = 5 * sig_
        rand = float(norm_val(0,max_,sig_))
        tot_err = (rand**2+scat**2)**.5

        return tot_err

    def t_(row):

        l_44 = row['log_L_agn_44']
        t_raw = tau_lum(l_44,m,c)
        return t_raw

    def pert(row,column):

        column = column[0]
        if column == 0: 
            raw = row['t_raw']
            err = row['log_t_c_err']

        elif column ==1:
            raw = row['log_L_agn_r']
            err = row['log_L_agn_err']


#       pertibation generation 
        min_ = -raw
        max_ =  raw

        pertibation = norm_val(min_,max_,err)
        val = raw + pertibation
        return val

    def lin_t(row):
        
        return 10**row['log_t_cent']

    def lin_t_err(row):

        err = 10**row['log_t_cent']*np.log(10.)*row['log_t_c_err']
        return err

    def frame():
        
        print 'Simulating the luminosity of',sim_n, 'objects:'
        draft = pd.DataFrame()
        draft_ = [[draft,draft]]

        draft['z'] = redshift_dist()
        func_app(dist,[['log_L_agn_r']],draft_)
        func_app(l_err,[['log_L_agn_err']],draft_)
        func_app(pert,[['log_L_agn',[1]]],draft_)
        print get_time(), 'to simulate every luminosity.'
        func_app(L_44,[['log_L_agn_44']],draft_)
        func_app(t_err,[['log_t_c_err']],draft_)
        func_app(t_,[['t_raw']],draft_)
        func_app(pert,[['log_t_cent',[0]]],draft_)
        func_app(lin_t,[['t_cent']],draft_)
        func_app(lin_t_err,[['t_c_err']],draft_)
        draft['t_c_nerr'] = draft['t_c_err']
        draft['t_c_perr'] = draft['t_c_err']
        draft = draft[(draft.log_t_cent < sim_c['t_max'])&(draft.log_t_cent > sim_c['t_min'])&(draft.log_L_agn > sim_c['l_min'])]
        print get_time(), 'to simulate everything.'
        print 'There are',len(draft.index),'points left inbounds'
        return draft


    return frame()


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
    """
    sim = simgen(sim_mode)
    L0 = sim['log_L_agn_44']
    t0 = sim['log_t_cent']
    L_err0 = sim['log_L_agn_err']
    t_err0 = sim['log_t_c_err'] 
    """
    L1 = full1['log_L_agn_44']
    t1 = full1['log_t_cent']
    L_err1 = full1['log_L_agn_err']
    t_err1 = [full1['log_t_c_nerr'],full1['log_t_c_perr']] 

    L2 = full2['log_L_agn_44']
    t2 = full2['log_t_cent']
    L_err2 = full2['log_L_agn_err']
    t_err2 = [full2['log_t_c_nerr'],full2['log_t_c_perr']] 

    L3 = full3['log_L_agn_44']
    t3 = full3['log_t_cent']
    L_err3 = full3['log_L_agn_err']
    t_err3 = [full3['log_t_c_nerr'],full3['log_t_c_perr']] 

    
    #plot labels
    Label0 = r'Log L vs Log $\tau$ Simulated Data'
    Label1 = r'Log L vs Log $\tau$ Bentz'
    Label2 = r'Log L vs Log $\tau$ SDSSRM'
    Label3 = r'Log L vs Log $\tau$ Du'
    Label4 = 'Best fit using ODR of all data'
    Label5 = 'Best fit using Bentz et. al. numbers'

    #pl.errorbar(L0, t0,xerr=L_err0, yerr = t_err0, fmt = 'y,', ecolor='y'  , elinewidth= .5, label=Label0,zorder=10)
    pl.errorbar(L1, t1,xerr=L_err1, yerr = t_err1, fmt = 'k.', ecolor='b'  , elinewidth= .5, label=Label1,zorder=10)
    pl.errorbar(L2, t2,xerr=L_err2, yerr = t_err2, fmt = 'b.', ecolor='k'  , elinewidth= .5, label=Label2,zorder=10)
    pl.errorbar(L3, t3,xerr=L_err3, yerr = t_err3, fmt = 'g.', ecolor='r'  , elinewidth= .5, label=Label3,zorder=10)
    pl.plot(x0,y0,'b:', label=Label4, linewidth=1.2)
    pl.plot(x1,y1,'g--',label=Label5, linewidth= .6)

    legend = pl.legend(loc="upper right")
    pl.set_xlabel(r'Log (Luminosity / $10^{44}$)', fontsize = 14.5)
    pl.set_ylabel('Log Tau', fontsize = 14.5)

    text = str('Regression slope is '+str(round(m,3))+ ' \ny-intercept is '+str(round(c,3)))
    plt.text(0.2,.5,text)

    pl.axis([-3,3,0,3])
    fig.savefig('%d_obsz_%d'%(sim_n,sim_mode),bbox_inches='tight',format='pdf')
    plt.show()

#L_tau_plot() 